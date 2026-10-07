"""
SAF-ThermoTwin  –  saf_thermotwin.py
=====================================
Validated Non-Ideal Brayton Cycle Framework for Thermodynamic Performance
and CORSIA-Aligned Lifecycle CO₂ Analysis of SAF Blends in Micro-Turbojet Engines.

Engine model : AMT Olympus HP (single-spool micro-turbojet)
Calibration  : JP-4 SLS test data, Leylek (2012) DSTO-TR-2757
Thermodynamics: Cantera (NASA polynomials, thermally perfect gas)
SAF pathways : HEFA-SPK, FT-SPK, ATJ-SPK, HFS-SIP  (ASTM D7566)

Workflow (adımlar / steps)
--------------------------
  Step 1  – Design-point data (AMT Olympus HP @ SLS)
  Step 2  – T4 calibration via bisection (target: EGT from Leylek 2012)
  Step 3  – Synthetic performance map calibration (Leylek operating line)
  Step 4  – SLS validation with synthetic maps
  Step 5a – Cruise baseline (10,000 m, Mach 0.6, JP-4)
  Step 5b – Cruise baseline (10,000 m, Mach 0.6, Jet-A)
  Step 6  – SAF comparison at maximum ASTM D7566 blend fractions

Not / Note
----------
Konsol çıktıları Türkçe yazılmıştır.
Console output is currently in Turkish.

Co-processing yakıtı ASTM D7566 kapsamında tanımlı bir yol değildir
ve makale sonuçlarına dahil edilmemiştir; yalnızca göstergeseldir.
The "Co-processing" fuel entry is NOT an ASTM D7566 pathway and is
NOT included in the paper results; its lifecycle value is indicative only.

Atıf / Citation
---------------
B. E. Emin, "SAF-ThermoTwin: Validated Non-Ideal Brayton Cycle Framework
for Thermodynamic Performance and Emissions Modelling of SAF Blends in
Micro Turbojet Engines," 2026 11th Int. Conf. on Recent Advances in Air
and Space Technologies (RAST), Istanbul, Turkiye, 2026, pp. 1–6.
DOI: https://doi.org/10.1109/RAST69551.2026.11672525

Lisans / License
----------------
Copyright (C) 2026  Berkay Eren Emin
This program is free software: you can redistribute it and/or modify
it under the terms of the GNU General Public License as published by
the Free Software Foundation, either version 3 of the License, or
(at your option) any later version.

This program is distributed in the hope that it will be useful,
but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
GNU General Public License for more details.

You should have received a copy of the GNU General Public License
along with this program. If not, see <https://www.gnu.org/licenses/>.

SPDX-License-Identifier: GPL-3.0-or-later
"""

import numpy as np
import cantera as ct
from scipy.interpolate import interp1d
from scipy.optimize import curve_fit

FUEL_DATABASE = {

    "JP-4": {
        "HC_ratio": 2.0,
        "LHV": 43.323,
        "density": 0.760,
        "CO2_factor": 3.140,
        "bio_carbon_fraction": 0.0,
        "lifecycle_CO2_gMJ": 89.0,
    },

    "Jet-A": {
        "formula_note": "Surrogate: C11H21 (n-decane/trimethylbenzene blend)",
        "HC_ratio": 1.909,
        "carbon_mass_fraction": 0.859,
        "hydrogen_mass_fraction": 0.136,
        "LHV": 43.2,                   # MJ/kg  – ASTM D7566-21 Tablo 1 [4]
        "density": 0.800,              # kg/L @ 15 °C – ASTM D7566-21 Tablo 1 [4]
        "CO2_factor": 3.149,
        "bio_carbon_fraction": 0.0,
        "lifecycle_CO2_gMJ": 89.0,     # gCO2e/MJ – ICAO CORSIA referans değeri [3]
        "description": "Conventional petroleum-derived jet fuel",
        "color": "#2c3e50",
    },
 
    "HEFA-SPK": {
        "formula_note": "Surrogate: C12H26 (iso-dodecane heavy)",
        "HC_ratio": 2.167,
        "carbon_mass_fraction": 0.849,
        "hydrogen_mass_fraction": 0.151,
        "LHV": 44.1,                   # MJ/kg  – ASTM D7566-21 Tablo A2.1 [4]
        "density": 0.758,              # kg/L @ 15 °C – ASTM D7566-21 Tablo A2.1
                                       #   (730–772 kg/m³ aralığı; [4])
        "CO2_factor": 3.113,
        "bio_carbon_fraction": 1.0,
        # ÖNCEKİ DEĞER: 15.0  →  yalnızca jatropha en iyi senaryosunu yansıtıyordu
        # YENİ DEĞER   : 22.5  →  CORSIA HEFA medyanı (ILUC dahil) [3]
        #   Aralık: 10.4 (jatropha) – 87.8 (palm) gCO2e/MJ
        #   Ortalama: 39.4 gCO2e/MJ, Medyan: 22.5 gCO2e/MJ
        "lifecycle_CO2_gMJ": 22.5,
        "description": "HEFA-SPK from bio-lipids (tallow, camelina, jatropha)",
        "lifecycle_note": (
            "CORSIA HEFA medyanı (ILUC dahil). "
            "Feedstock'a göre 10.4–87.8 gCO2e/MJ arasında değişir. "
            "Kaynak: Kurzawska-Pietrowicz (2023) [3]."
        ),
        "color": "#27ae60",
    },
 
    "FT-SPK": {
        "formula_note": "Surrogate: C10H22 (n-decane)",
        "HC_ratio": 2.200,
        "carbon_mass_fraction": 0.843,
        "hydrogen_mass_fraction": 0.157,
        "LHV": 44.3,                   # MJ/kg  – Rye et al. (2010) [1]; ASTM D7566-21 [4]
        "density": 0.740,              # kg/L @ 15 °C – ASTM D7566-21 Tablo A1.1
                                       #   (730–770 kg/m³ aralığı; [4])
        "CO2_factor": 3.092,
        "bio_carbon_fraction": 1.0,
        # DEĞER KORUNDU: 10 gCO2e/MJ
        #   FT-SPK medyanı: 6.8 gCO2e/MJ, ortalaması: 3.5 gCO2e/MJ [3]
        #   Atık biyokütle senaryoları: 7.7–8.3 gCO2e/MJ [3]
        #   10 gCO2e/MJ bu aralık için makul bir üst tahmin değeridir.
        "lifecycle_CO2_gMJ": 10.0,
        "description": "FT-SPK from biomass gasification",
        "lifecycle_note": (
            "CORSIA FT-SPK aralığı: -22.5 ile +20.8 gCO2e/MJ. "
            "Medyan: 6.8 gCO2e/MJ. 10 değeri orta-konservatif senaryo. "
            "Kaynak: Kurzawska-Pietrowicz (2023) [3]; Hileman & Stratton (2014) [2]."
        ),
        "color": "#2980b9",
    },
 
    "ATJ-SPK": {
        "formula_note": "Surrogate: C9H20 (iso-nonane)",
        "HC_ratio": 2.222,
        "carbon_mass_fraction": 0.838,
        "hydrogen_mass_fraction": 0.162,
        "LHV": 43.9,                   # MJ/kg  – ASTM D7566-21 Tablo A5.1 [4]
        "density": 0.755,              # kg/L @ 15 °C – ASTM D7566-21 Tablo A5.1
                                       #   (730–770 kg/m³ aralığı; [4])
        "CO2_factor": 3.074,
        "bio_carbon_fraction": 1.0,
        # ÖNCEKİ DEĞER: 38.0  →  literatürün biraz üzerindeydi
        # YENİ DEĞER   : 33.0  →  izobutanol ATJ medyanı (31.3) ile
        #                          etanol ATJ medyanı (32.5) ortalaması [3]
        #   İzobutanol aralığı: -10.7 – +85.5 gCO2e/MJ, medyan: 31.3 [3]
        #   Etanol aralığı    : -6.8  – +100.6 gCO2e/MJ, medyan: 32.5 [3]
        "lifecycle_CO2_gMJ": 33.0,
        "description": "ATJ-SPK from alcohol dehydration",
        "lifecycle_note": (
            "İzobutanol ATJ medyanı: 31.3 gCO2e/MJ; etanol ATJ medyanı: 32.5 gCO2e/MJ. "
            "33.0 her iki yolun ağırlıklı ortalamasıdır. "
            "Kaynak: Kurzawska-Pietrowicz (2023) [3]."
        ),
        "color": "#8e44ad",
    },
 
    "HFS-SIP": {
        "formula_note": "Surrogate: C15H32 (farnesane)",
        "HC_ratio": 2.133,
        "carbon_mass_fraction": 0.851,
        "hydrogen_mass_fraction": 0.149,
        "LHV": 43.6,                   # MJ/kg  – ASTM D7566-21 Tablo A3.1 (min 43.5) [4]
        # ASTM D7566-21 Tablo A3.1: 765–780 kg/m³ → 0.770 kg/L [4]
        # ÖNCEKİ DEĞER: 0.763  →  spec. minimumunun (765 kg/m³) altındaydı
        "density": 0.770,              # kg/L @ 15 °C – ASTM D7566-21 Tablo A3.1 orta değer [4]
        "CO2_factor": 3.120,
        "bio_carbon_fraction": 1.0,
        "max_blend_fraction": 0.10,    # %10 karışım sınırı – ASTM D7566-21 Madde 6.1.3 [4]
        # ÖNCEKİ DEĞER: 25.0  →  literatürün önemli ölçüde altındaydı
        # YENİ DEĞER   : 44.0  →  CORSIA SIP medyanı [3]
        #   Şeker kamışı (Brezilya): ~44.1 gCO2e/MJ
        #   Şeker pancarı (AB)     : ~52.6 gCO2e/MJ
        #   Medyan: 44 gCO2e/MJ, Ortalama: 46 gCO2e/MJ
        "lifecycle_CO2_gMJ": 44.0,
        "description": "HFS-SIP from fermented sugars (farnesane, max 10%)",
        "lifecycle_note": (
            "CORSIA SIP medyanı: 44 gCO2e/MJ (ILUC dahil). "
            "Şeker kamışı (BR): ~44.1, şeker pancarı (AB): ~52.6 gCO2e/MJ. "
            "Kaynak: Kurzawska-Pietrowicz (2023) [3]."
        ),
        "color": "#e67e22",
    },
 
    "Co-processing": {
        "formula_note": "Intermediate between Jet-A and HEFA",
        "HC_ratio": 2.000,
        "carbon_mass_fraction": 0.854,
        "hydrogen_mass_fraction": 0.143,
        "LHV": 43.5,
        "density": 0.780,
        "CO2_factor": 3.130,
        "bio_carbon_fraction": 0.30,
        # Co-processing, ASTM D7566-21 kapsamında tanımlı bir yol değildir.
        # 55 gCO2e/MJ değeri %30 biyojenik karbon fraksiyonuna orantılı tahmindir:
        #   0.70 × 89 + 0.30 × 22.5 ≈ 62 + 6.75 ≈ 69 gCO2e/MJ
        # Not: D7566-21 kapsamında onaylı bir yol olmadığından bu değer
        # yalnızca göstergesel kabul edilmelidir. [1,2,4]
        "lifecycle_CO2_gMJ": 55.0,
        "description": "Co-processed bio-lipid + petroleum fraction (D7566 dışı yol)",
        "lifecycle_note": (
            "ASTM D7566-21 kapsamında tanımlı değildir. "
            "55 gCO2e/MJ değeri orantılı karışım tahminidir; "
            "kesin LCA için bağımsız doğrulama önerilir. "
            "Kaynak: Rye et al. (2010) [1]; Hileman & Stratton (2014) [2]."
        ),
        "color": "#c0392b",
    },
}

# =============================================================================
# KARIŞIM HESAPLAYICI
# =============================================================================
def compute_blend_properties_jp4(saf_type, blend_fraction, verbose=True):
    #jeta = FUEL_DATABASE["Jet-A"]
    jeta = FUEL_DATABASE["JP-4"]
    saf  = FUEL_DATABASE[saf_type]

    max_blend = saf.get("max_blend_fraction", 0.50)
    if blend_fraction > max_blend:
        print(f"  [UYARI] {saf_type} ASTM karışım sınırı {max_blend*100:.0f}%. "
              f"İstenen {blend_fraction*100:.0f}% sınırı aşıyor.")

    rho_jeta = jeta["density"]
    rho_saf  = saf["density"]
    mass_jeta = (1 - blend_fraction) * rho_jeta
    mass_saf  = blend_fraction * rho_saf
    total_mass = mass_jeta + mass_saf
    x_jeta = mass_jeta / total_mass
    x_saf  = mass_saf  / total_mass

    blend = {
        "saf_type":            saf_type,
        "vol_fraction_saf":    blend_fraction,
        "mass_fraction_saf":   x_saf,
        "LHV":                 x_jeta * jeta["LHV"] + x_saf * saf["LHV"],
        "density":             1.0 / (x_jeta / rho_jeta + x_saf / rho_saf),
        "HC_ratio":            x_jeta * jeta["HC_ratio"] + x_saf * saf["HC_ratio"],
        "CO2_factor":          x_jeta * jeta["CO2_factor"] + x_saf * saf["CO2_factor"],
        "bio_carbon_fraction": x_jeta * jeta["bio_carbon_fraction"] + x_saf * saf["bio_carbon_fraction"],
        "lifecycle_CO2_gMJ":   x_jeta * jeta["lifecycle_CO2_gMJ"] + x_saf * saf["lifecycle_CO2_gMJ"],
    }
    blend["net_CO2_factor"] = blend["CO2_factor"] * (1 - blend["bio_carbon_fraction"])

    if verbose:
        print(f"\n  Karışım: {blend_fraction*100:.1f}% {saf_type} + {(1-blend_fraction)*100:.1f}% Jet-A (hacimce)")
        print(f"    LHV                 : {blend['LHV']:.3f} MJ/kg")
        print(f"    Yoğunluk            : {blend['density']:.4f} kg/L")
        print(f"    H/C oranı           : {blend['HC_ratio']:.4f}")
        print(f"    Stoik. CO₂ faktörü  : {blend['CO2_factor']:.4f} kg CO₂/kg yakıt")
        print(f"    Biyojenik C fraks.  : {blend['bio_carbon_fraction']*100:.1f}%")
        print(f"    Net fosil CO₂       : {blend['net_CO2_factor']:.4f} kg CO₂/kg yakıt")
        print(f"    Yaşam döngüsü CO₂   : {blend['lifecycle_CO2_gMJ']:.1f} gCO₂eq/MJ")
    return blend

def compute_blend_properties(saf_type, blend_fraction, verbose=True):
    jeta = FUEL_DATABASE["Jet-A"]
    saf  = FUEL_DATABASE[saf_type]

    max_blend = saf.get("max_blend_fraction", 0.50)
    if blend_fraction > max_blend:
        print(f"  [UYARI] {saf_type} ASTM karışım sınırı {max_blend*100:.0f}%. "
              f"İstenen {blend_fraction*100:.0f}% sınırı aşıyor.")

    rho_jeta = jeta["density"]
    rho_saf  = saf["density"]
    mass_jeta = (1 - blend_fraction) * rho_jeta
    mass_saf  = blend_fraction * rho_saf
    total_mass = mass_jeta + mass_saf
    x_jeta = mass_jeta / total_mass
    x_saf  = mass_saf  / total_mass

    blend = {
        "saf_type":            saf_type,
        "vol_fraction_saf":    blend_fraction,
        "mass_fraction_saf":   x_saf,
        "LHV":                 x_jeta * jeta["LHV"] + x_saf * saf["LHV"],
        "density":             1.0 / (x_jeta / rho_jeta + x_saf / rho_saf),
        "HC_ratio":            x_jeta * jeta["HC_ratio"] + x_saf * saf["HC_ratio"],
        "CO2_factor":          x_jeta * jeta["CO2_factor"] + x_saf * saf["CO2_factor"],
        "bio_carbon_fraction": x_jeta * jeta["bio_carbon_fraction"] + x_saf * saf["bio_carbon_fraction"],
        "lifecycle_CO2_gMJ":   x_jeta * jeta["lifecycle_CO2_gMJ"] + x_saf * saf["lifecycle_CO2_gMJ"],
    }
    blend["net_CO2_factor"] = blend["CO2_factor"] * (1 - blend["bio_carbon_fraction"])

    if verbose:
        print(f"\n  Karışım: {blend_fraction*100:.1f}% {saf_type} + {(1-blend_fraction)*100:.1f}% Jet-A (hacimce)")
        print(f"    LHV                 : {blend['LHV']:.3f} MJ/kg")
        print(f"    Yoğunluk            : {blend['density']:.4f} kg/L")
        print(f"    H/C oranı           : {blend['HC_ratio']:.4f}")
        print(f"    Stoik. CO₂ faktörü  : {blend['CO2_factor']:.4f} kg CO₂/kg yakıt")
        print(f"    Biyojenik C fraks.  : {blend['bio_carbon_fraction']*100:.1f}%")
        print(f"    Net fosil CO₂       : {blend['net_CO2_factor']:.4f} kg CO₂/kg yakıt")
        print(f"    Yaşam döngüsü CO₂   : {blend['lifecycle_CO2_gMJ']:.1f} gCO₂eq/MJ")
    return blend

# =============================================================================
# İDEAL OLMAYAN BRAYTON ÇEVRİMİ – TURBOJET
# =============================================================================
def actual_brayton_cycle_turbojet(
    h, M_ambient, rp_comp, eta_comp, eta_turb, eta_nozz, T4,
    fuel_blend, deltaP_comb=0.05, eta_mech=0.99, eta_comb=0.98,
    mechanism="air.yaml",
):
    """
    Cantera tabanlı gerçek hava termodinamiği ile ideal olmayan Brayton çevrimi.
    Yakıt termokimyası: fuel_blend dict üzerinden.

    Parametreler
    ------------
    h           : Uçuş irtifası [m]
    M_ambient   : Uçuş Mach sayısı [-]
    rp_comp     : Kompresör basınç oranı [-]
    eta_comp    : İzentropik kompresör verimi [-]
    eta_turb    : İzentropik türbin verimi [-]
    eta_nozz    : Nozul verimi [-]
    T4          : Türbin giriş sıcaklığı [K]
    fuel        : Yakıt sözlüğü (JETA_FUEL/SAF)
    deltaP_comb : Yanma odası basınç kaybı (fraksiyonel) [-]
    eta_mech    : Mekanik verim [-]
    eta_comb    : Yanma verimi [-]
    """

    import warnings

    gas = ct.Solution(mechanism)

    def state_from_gas(g):
        return {
            "T":     g.T,
            "P":     g.P / 1e3,
            "h":     g.h,
            "s":     g.s,
            "v":     1.0 / g.density,
            "cp":    g.cp,
            "gamma": g.cp / g.cv,
        }

    # ── ISA Atmosfer ─────────────────────────────────────────────────────────
    def ISA_atmosphere(altitude):
        g = 9.80665
        gas.TP = 288.15, 101325.0
        R = ct.gas_constant / gas.mean_molecular_weight
        T_sl, P_sl = 288.15, 101325.0

        # Atmosferin deniz seviyesinden 11 km'ye kadar olan kısmına Troposfer denir. 
        #   Tüm hava olayları ve ticari uçuşların büyük kısmı burada gerçekleşir. 
        #   Troposferde yukarı çıktıkça sıcaklık her 1 kilometrede 6.5 C düşer. 
        #   Sıcaklık (T) buna göre doğrusal olarak hesaplanır. Sıcaklık doğrusal olarak 
        #   düştüğü için basınç hidrostatik denge denklemine göre polytropik bir eğri 
        #   izleyerek düşer. Bu karmaşık üslü denklem, akışkanlar mekaniğindeki barometrik formülün ta kendisidir.
        if altitude <= 11000:
            L = -0.0065
            T = T_sl + L * altitude
            P = P_sl * (T / T_sl) ** (-g / (L * R))

        # 11 km irtifaya ulaştığımızda atmosferin karakteri aniden değişir. Burası 
        #   Tropopoz (Tropopause) sınırıdır ve Stratosfer başlar. 11 km ile 20 km arasında 
        #   sıcaklık düşmeyi bırakır ve yaklaşık -56.5 C (216.65 K) seviyesinde sabit kalır. 
        #   Sıcaklık artık sabit (izotermal) olduğu için, barometrik denklem değişir ve 
        #   basınç artık üslü değil, üstel (eksponansiyel, e^-z) olarak düşmeye başlar.    
        elif altitude <= 20000:
            # Stratosfer: sabit sıcaklık katmanı (ISA, 11–20 km)
            L_trop = -0.0065  # Troposfer ısıl derecesi [K/m] – yalnızca troposfer sınırı için
            T11 = T_sl + L_trop * 11000  # = 216.65 K
            P11 = P_sl * (T11 / T_sl) ** (-g / (L_trop * R))
            T = T11
            P = P11 * np.exp(-g * (altitude - 11000) / (R * T))
        else:
            raise ValueError("ISA modeli yalnızca 20 km'ye kadar geçerlidir.")
        gas.TP = T, P
        return {"T_ambient": T, "P_ambient": P, "rho_ambient": gas.density}

    # ── Durum 1: Ortam / Difüzör Girişi ─────────────────────────────────────
    amb = ISA_atmosphere(h)
    T_static, P_static = amb["T_ambient"], amb["P_ambient"]
    gas.TP = T_static, P_static
    state1 = state_from_gas(gas)
    a0 = gas.sound_speed
    V0 = M_ambient * a0

    # ── Durum 2: Difüzör Çıkışı ─────────────────────────────────────────────
    h0 = state1["h"] + 0.5 * V0 ** 2
    gamma = gas.cp / gas.cv

    # Havayı motora göre durma noktasına getirdiğimizde, sahip olduğu statik 
    #   entalpi (h) ile kinetik enerjisi (V^2/2) toplanarak "toplam (durgunluk) 
    #   entalpisi" (h0) elde edilir. T0, bu toplam enerjinin sıcaklık karşılığıdır.
    #   P0_ideal ise bu yavaşlatma ve sıkıştırma işleminin %100 verimle (hiçbir 
    #   kayıp olmadan, izantropik) gerçekleşmesi durumunda ulaşılacak teorik 
    #   maksimum basınçtır.
    T0 = T_static * (1 + (gamma - 1) / 2 * M_ambient ** 2)
    P0_ideal = P_static * (T0 / T_static) ** (gamma / (gamma - 1))

    # Yüksek hızlarda ilerleyen bir hava aracında, hava motorun girişine adeta "çarpıp" sıkışır. 
    #   Herhangi bir hareketli parça (kompresör) kullanmadan, sadece geometrik daralma ve hızın 
    #   düşürülmesiyle elde edilen bu doğal basınç artışına Ram Etkisi denir. Ramjet motorlarının 
    #   çalışma prensibi tamamen buna dayanır; turbojetlerde ise yüksek hızlarda kompresörün 
    #   işini inanılmaz derecede kolaylaştırır.

    # Özetle Ram etkisi; aracın hareketinden gelen kinetik enerjinin, akış yavaşlatılarak basınca çevrilmesi sonucu oluşan doğal sıkıştırmadır.
    if M_ambient <= 1.0:

        # Şok dalgası yoktur. Kayıplar sadece difüzör yüzeyindeki vizkoz sürtünmelerden 
        #   kaynaklanır. Yüzde 2'lik bir toplam basınç kaybı makul bir yaklaşımdır.

        # Hava alığı (intake/diffuser) kayıpları, toplam basınç geri kazanım katsayısı (ram recovery factor) 
        #   kullanılarak modellenmiştir. Bu yaklaşımda difüzör çıkışındaki toplam basınç, serbest akış toplam 
        #   basıncının bir fraksiyonu olarak ifade edilmiştir.
        #   Literatürde, özellikle süpersonik akışlarda şok dalgalarından kaynaklanan toplam basınç kayıplarını 
        #   temsil etmek amacıyla bu yaklaşım yaygın olarak kullanılmaktadır. Bu çalışmada da süpersonik rejim 
        #   için MIL-spec tabanlı ampirik bir ram recovery bağıntısı kullanılmıştır.
        #   Subsonik rejim için literatürde alternatif olarak izentropik verim (η_intake) tabanlı bir yaklaşım 
        #   da mevcuttur. Bu yöntemde toplam basınç artışı, izentropik sıkıştırma ilişkisine bir verim katsayısı 
        #   eklenerek hesaplanmaktadır (örneğin kitapta verilen Eq. 3.10b). Söz konusu kaynakta ayrıca ram recovery 
        #   ile izentropik verim yaklaşımlarının subsonik rejimde approximate olarak birbirine eşdeğer olduğu 
        #   (η_r ≈ η_i) ifade edilmektedir.
        #   Bu nedenle, bu çalışmada tüm hız aralığında tek ve tutarlı bir model kullanmak amacıyla ram recovery 
        #   yaklaşımı tercih edilmiştir. Subsonik rejimde iki yöntem arasındaki farkın küçük olduğu (≈%1 mertebesinde)
        #   bilindiğinden, bu seçimin sonuçlar üzerindeki etkisinin ihmal edilebilir olduğu kabul edilmiştir.
        ram_recovery = 0.98
    elif M_ambient <= 5.0:

        # Hava sesten hızlıdır. Motorun havayı yutabilmesi için yavaşlatması gerekir, 
        #   bu da difüzör içinde şok dalgaları (oblique and normal shocks) oluşturur. Şok dalgaları 
        #   ciddi entropi artışına ve dolayısıyla toplam basınç kaybına yol açar. Hız arttıkça şoklar 
        #   şiddetlenir ve geri kazanım (recovery) bu ampirik formüle göre logaritmik olarak düşer.
        #   (MIL-E-5007D — normal şok + oblique şok kombinasyonu)
        ram_recovery = 1.0 - 0.075 * (M_ambient - 1) ** 1.35
    else:

        # Mach 5 ve sonrasında şok dalgaları o kadar şiddetlidir ki, standart eğriler yetersiz kalır. 
        #   Bu bölgede toplam basınç kaybı devasa boyutlara ulaşır. Bu formül, hipersonik rejimdeki 
        #   bu dramatik düşüşü modelleyen standart bir yaklaşımdır.
        #   (NASA hipersonik) 
        ram_recovery = 800 / (M_ambient ** 4 + 935)

    P_diff_out = P0_ideal * ram_recovery
    gas.HP = h0, P_diff_out
    state2 = state_from_gas(gas)

    # ── Durum 3: Kompresör Çıkışı ────────────────────────────────────────────
    T2, P2 = state2["T"], state2["P"]
    P3_Pa = rp_comp * P2 * 1e3
    gas.TP = T2, P2 * 1e3
    h2, s2 = gas.h, gas.s

    gas.SP = s2, P3_Pa
    h3s = gas.h
    h3 = h2 + (h3s - h2) / eta_comp
    gas.HP = h3, P3_Pa
    state3 = state_from_gas(gas)
    w_comp = h3 - h2

    # ── Durum 4: Yanma Odası Çıkışı ──────────────────────────────────────────
    P4_Pa = (1.0 - deltaP_comb) * state3["P"] * 1e3
    gas.TP = state3["T"], P4_Pa
    h3_at_P4 = gas.h

    gas.TP = T4, P4_Pa
    state4 = state_from_gas(gas)
    h4 = state4["h"]
    
    LHV_Jkg = fuel_blend["LHV"] * 1e6
    q_in_air = h4 - h3_at_P4
    f = q_in_air / (eta_comb * LHV_Jkg)

    gas.TP = state4["T"], P4_Pa
    s4 = gas.s

    q_in = f * LHV_Jkg   # [J / kg_hava]

    # ── Durum 5: Türbin Çıkışı ────────────────────────────────────────────────
    
    # Turbojetlerde türbinin ana görevi kompresörü çevirmektir. Sisteme yanma 
    #   odasında yakıt eklendiği için türbinden geçen debi (1+f) oranında artar. 
    #   Gerekli türbin işi (w_turb_req), kompresör işi, kütle debisi artışı ve 
    #   şaft mekanik verimi (eta_mech) kullanılarak hesaplanır.
    
    #   Güç dengesi:  w_comp [J/kg_hava] = (1+f) · w_turb [J/kg_ürün] · η_mech
    w_turb_required = w_comp / ((1.0 + f) * eta_mech)

    # Burada kullanılan yöntem, İkiye Bölme (Bisection) algoritmasına dayanmaktadır. 
    #   Termodinamik özellikler (özellikle Cantera gibi gerçek gaz kütüphanelerinde) doğrusal 
    #   değildir. "Bana tam olarak "gereken türbin işi" kadar iş verecek çıkış basıncı P_5 nedir?" 
    #   sorusunun doğrudan, tek satırlık bir cebirsel cevabı yoktur. 

    #   Bu yüzden kod bir bisection (ikiye bölme) algoritması kuruyor: 
    #   - Mantıklı bir alt (P5_low) ve üst (P5_high) basınç sınırı belirliyor.
    #   - İkisinin tam ortasını (P5_trial) deniyor.
    #   - Eğer bu deneme basıncında üretilen iş, gereken işten azsa (basınç yeterince düşmemiş demektir), 
    #       üst sınırı aşağı çekiyor. Fazlaysa, alt sınırı yukarı çekiyor.
    #   - İstenen işe çok yaklaşana kadar (abs < 0.5) veya 60 iterasyon dolana kadar 
    #       bu deneme-yanılmayı sistematik olarak sürdürüyor.

    # Türbinin amacı kompresörü çevirecek işi üretmektir. Termodinamik olarak bir gazdan 
    # mekanik iş çekebilmek için gazın genişlemesi, yani basıncının düşmesi şarttır. (P5_high = P4_Pa * 0.99) 

    # Türbinden çıkan gaz (Durum 5), lüleye (nozzle) girer.  Lülenin görevi, bu sıcak gazı dışarıdaki 
    # ortam (atmosfer) basıncına ($P_{static}$) kadar genişleterek hızlandırmak ve uçağı itecek o devasa 
    # tepki kuvvetini (thrust) üretmektir. Lülenin gazı dışarı atabilmesi için, türbinden çıkan gazın 
    # basıncının (P_5), dışarıdaki atmosfer basıncından daha yüksek olması gerekir. (P5_low  = P_static * 0.5)
    
    P5_low  = P_static * 0.5
    P5_high = P4_Pa * 0.99
    w_turb_actual = 0.0
    
    _bisect_converged = False
    for _ in range(60):
        P5_trial = (P5_low + P5_high) / 2
        gas.SP = s4, P5_trial
        h5s = gas.h
        w_turb_actual = eta_turb * (h4 - h5s)
        if abs(w_turb_actual - w_turb_required) < 0.5:  # J/kg toleransı
            P5_Pa = P5_trial
            _bisect_converged = True
            break
        if w_turb_actual < w_turb_required:
            P5_high = P5_trial
        else:
            P5_low = P5_trial
    else:
        P5_Pa = P5_trial

        # Eğer kod buraya düşerse, matematiksel bir sıkıntı var demektir 
        #   (örneğin başlangıç sınırları hatalı seçilmiş olabilir veya o şartlarda o işi 
        #   üretmek fiziksel olarak imkansız olabilir). Bu durumda kod çökmek yerine son 
        #   denediği değeri kabul eder ama kullanıcıya (örneğin bir termodinamik motor 
        #   modeli simülasyonunda çok kritik olan) "Dikkat, yakınsamadı, hesaplarda 
        #   şu kadarlık hata var" uyarısı verir.
        warnings.warn(
            f"[SAF-ThermoTwin] Türbin P5 biseksiyon yakınsamadı! "
            f"Artık hata = {abs(w_turb_actual - w_turb_required):.1f} J/kg. "
            f"Giriş: PR={rp_comp}, T4={T4} K, eta_turb={eta_turb:.3f}",
            RuntimeWarning, stacklevel=3
        )

    # Gerçek türbin çıkış durumu (irreversibility dahil)
    gas.SP = s4, P5_Pa
    h5s = gas.h
    h5 = h4 - eta_turb * (h4 - h5s)

    gas.HP = h5, P5_Pa
    state5 = state_from_gas(gas)
    rp_turb = P4_Pa / P5_Pa
    w_turb = h4 - h5

    # ── Durum 6: Nozul Çıkışı ─────────────────────────────────────────────────
    # Ortam basıncı (nozulun genişlediği hedef basınç)
    P_ambient = state1["P"] * 1e3

    # Türbin çıkışındaki durum (nozul girişi)
    gas.TP = state5["T"], P5_Pa
    h5_gas, s5 = gas.h, gas.s

    # Gazın anlık özgül ısı oranı (gamma)
    gamma5 = gas.cp / gas.cv

    # Öncelikle türbin çıkışındaki gazın özgül ısı oranı (gamma) hesaplanıyor. 
    #   Ardından, sıkıştırılabilir akış (compressible flow) teorisinin en meşhur 
    #   denklemlerinden biri kullanılarak Kritik Basınç Oranı (crit) bulunuyor.

    # Bir gazı, daralan bir kanalda (converging nozzle) hızlandırdığında 
    #   ulaşabileceği maksimum hız ses hızıdır (Mach 1). Akışın tam çıkış noktasında 
    #   (boğazda) Mach 1'e ulaşabilmesi için, içerideki basıncın (P_5), 
    #   dışarıdaki basınca (P_ambient) oranının bu crit değerinden büyük veya ona eşit olması gerekir.
    crit = 1.0 / (1.0 - (1.0/eta_nozz) * ((gamma5 - 1.0)/(gamma5 + 1.0))) \
       ** (gamma5 / (gamma5 - 1.0))
    choked = (P5_Pa / P_ambient) >= crit #Turbine exit flow velocity is assumed negligible, therefore total and static pressures are considered equal P5_Pa = P5_total_Pa = P5_static_Pa
    # Türbin çıkışındaki akış hızının, toplam enerji dengesi içerisindeki katkısının ihmal edilebilir olduğu varsayılmıştır. Bu kabul doğrultusunda türbin çıkışındaki toplam basınç (stagnation pressure) ile statik basınç eşit kabul edilmiştir.
    #   Fiziksel olarak, türbin içerisindeki akışın temel amacı, gazın genleşmesi yoluyla mekanik iş üretmektir. Bu süreçte akışın büyük bir kısmı, kinetik enerjiye dönüşmek yerine türbin kanatları üzerinde iş üretmek için kullanılır. Dolayısıyla türbin çıkışındaki hızlar, nozul çıkışına kıyasla oldukça düşüktür.
    #   Sıkıştırılabilir akış teorisine göre toplam basınç ile statik basınç arasındaki fark, akışın Mach sayısına bağlıdır. Düşük Mach sayılarında bu fark oldukça küçük olduğundan, türbin çıkışında toplam ve statik büyüklüklerin eşit kabul edilmesi yaygın bir mühendislik yaklaşımıdır.
    #   Bu varsayım, gaz türbini çevrim analizlerinde sıklıkla kullanılmakta olup, nozul giriş koşullarının belirlenmesinde hesaplamayı basitleştirirken sonuçlar üzerinde ihmal edilebilir bir hata oluşturmaktadır.

    # Eğer elimizdeki basınç oranı, kritik orandan büyükse veya ona eşitse, 
    #   akış "choked" yani boğulmuş durumdadır. Bu durumda gaz, nozulun dar 
    #   noktasında Mach 1'e ulaşır ve daha fazla hızlanamaz. Nozulun çıkışında 
    #   basınç, kritik basınç oranına göre belirlenir. Eğer choked ise, P6_Pa 
    #   kritik basınç oranına göre hesaplanır; değilse, gaz tamamen 
    #   genişleyerek dış ortam basıncına eşit olur.

    # ─────────────────────────────────────────────────────────────────────────
    # NOZUL ÇIKIŞ BASINCI
    #
    # Choked durumda:
    #   Gaz ambient basınca kadar genişleyemez
    #   Çıkış basıncı kritik basınç olur
    #
    # Unchoked durumda:
    #   Gaz tamamen genişler → P6 = P_ambient
    # ─────────────────────────────────────────────────────────────────────────

    if choked:
        P6_Pa = P5_Pa / crit
    else:
        P6_Pa = P_ambient


    # ─────────────────────────────────────────────────────────────────────────
    # İZENTROPİK GENİŞLEME (ideal nozzle)
    #
    # s = sabit alınarak (isentropic expansion)
    # h6s bulunur → bu "ideal" durum
    # ─────────────────────────────────────────────────────────────────────────
    gas.SP = s5, P6_Pa
    h6s = gas.h

    # ─────────────────────────────────────────────────────────────────────────
    # GERÇEK NOZUL (verim dahil)
    #
    # Nozul verimi:
    #   η_nozz = (gerçek kinetik enerji) / (ideal kinetik enerji)
    #
    # Entalpi formu:
    #   h6 = h5 - η_nozz * (h5 - h6s)
    #
    # Bu, irreversibility (kayıplar) içerir
    # ─────────────────────────────────────────────────────────────────────────
    h6 = h5_gas - eta_nozz * (h5_gas - h6s)
    gas.HP = h6, P6_Pa
    state6 = state_from_gas(gas)

    # ─────────────────────────────────────────────────────────────────────────
    # ÇIKIŞ HIZI (EN KRİTİK NOKTA)
    #
    # HER ZAMAN enerji denkleminden hesaplanır:
    #
    #   V = sqrt(2 * (h5 - h6))
    #
    # Neden?
    # Çünkü bu:
    #   → tüm termodinamik etkileri içerir
    #   → değişken cp, gerçek gaz vs. dahil
    #
    # NOT:
    # Mach=1 → V = sqrt(gamma*R*T) yaklaşımı burada kullanılmaz
    # çünkü bu sadece ideal ve sabit cp varsayımıdır
    # ─────────────────────────────────────────────────────────────────────────

    V6 = np.sqrt(max(0.0, 2.0 * (h5_gas - h6)))

    # ── Özgül Performans Metrikleri ───────────────────────────────────────────
    w_net = w_turb - w_comp
    eff_cycle = w_net / q_in if q_in > 0 else 0
    bwr = w_comp / ((1.0 + f) * w_turb)

    rho6 = gas.density  
    if abs(P6_Pa - P_ambient) > 1e-3 and V6 > 0.0:
        pressure_thrust = (P6_Pa - P_ambient) * (1.0 + f) / (rho6 * V6)
    else:
        pressure_thrust = 0.0
        
    # Özgül İtki (Specific Thrust)
    specific_thrust = (1.0 + f) * V6 - V0 + pressure_thrust

    SFC = (f / specific_thrust * 1e6) if specific_thrust > 0 else float("nan")

    # Termal verim: kinetik enerji artışı / eklenen ısı
    # Çıkış kinetik enerjisi: (1+f)·V6²/2 (yakıt kütlesi dahil)
    thermal_eff   = (0.5 * ((1.0 + f) * V6 ** 2 - V0 ** 2)) / q_in if q_in > 0 else 0
    # İtici verim: η_p = F·V0 / (ΔKE) = specific_thrust·V0 / (0.5·((1+f)·V6²−V0²))
    _dKE = 0.5 * ((1.0 + f) * V6 ** 2 - V0 ** 2)
    propulsive_eff = (specific_thrust * V0 / _dKE) if (_dKE > 0 and V0 > 0) else 0.0
    overall_eff   = thermal_eff * propulsive_eff

    CO2_per_kgair     = f * fuel_blend["CO2_factor"]
    net_CO2_per_kgair = f * fuel_blend["net_CO2_factor"]
    CO2_per_thrust_g_N     = (CO2_per_kgair * 1000 / specific_thrust) if specific_thrust > 0 else float("nan")
    net_CO2_per_thrust_g_N = (net_CO2_per_kgair * 1000 / specific_thrust) if specific_thrust > 0 else float("nan")

    return {
        "state1": state1, "state2": state2, "state3": state3,
        "state4": state4, "state5": state5, "state6": state6,
        "q_in": q_in, "q_out": state6["h"] - state1["h"],
        "w_comp": w_comp, "w_turb": w_turb, "w_net": w_net,
        "efficiency_cycle": eff_cycle, "bwr": bwr, "rp_turb": rp_turb,
        "V0": V0, "V6": V6, "choked": choked,
        "specific_thrust": specific_thrust,
        "SFC_mg_Ns": SFC, "FAR": f,
        "thermal_efficiency": thermal_eff,
        "propulsive_efficiency": propulsive_eff,
        "overall_efficiency": overall_eff,
        "CO2_per_kgair": CO2_per_kgair,
        "net_CO2_per_kgair": net_CO2_per_kgair,
        "CO2_per_thrust_g_N": CO2_per_thrust_g_N,
        "net_CO2_per_thrust_g_N": net_CO2_per_thrust_g_N,
        "fuel_blend": fuel_blend,
        # Verim bilgileri (harita veya sabit)
        "eta_comp_used": eta_comp,
        "eta_turb_used": eta_turb,
    }


# =============================================================================
# HİBRİT SENTETİK PERFORMANS HARİTASI
# =============================================================================

def calibrate_hybrid_map_from_leylek(
    PR_design=3.3812, eta_c_design=0.744
):
    """
    Leylek (2012) Fig.16, 18, 19'dan sayısallaştırılan operating line verisini
    kullanarak PR ve η_c harita parametrelerini belirler.
 
    Yöntem
    ------
    1. Corrected mass flow – corrected speed eşleşmesi (Fig.16)
    2. PR ve η_c'yi corrected mass flow üzerinden corrected speed'e bağla (Fig.19)
    3. Nr = N_corr / N_design olarak normalize et
    4. curve_fit: PR  = 1 + (PR_des-1)·Nr^a  → a belirlenir
                  η_c = η_c,des · exp(-kc·(Nr-1)²)  → kc belirlenir
    5. Fit kalitesi (R²) raporlanır
 
    Dönüş
    ------
    dict: a_fit, kc_fit, R2_PR, R2_eta, Nr_grid, PR_grid, eta_grid
    """
    # ── Sayısallaştırılmış Leylek Verileri ──────────────────────────────────
    # Corrected mass flow [kg/s] vs corrected speed [RPM]  (Fig.16)
    N_design_corr = 106064.0   # RPM – Fig.18 tepesinden
 
    mfr_raw = np.array([
        (50880,  0.21977), (62902,  0.27477), (68228,  0.30373),
        (73630,  0.33312), (77511,  0.35536), (81543,  0.37887),
        (84587,  0.39691), (87022,  0.40909), (89533,  0.42420),
        (91891,  0.43721), (93641,  0.44477), (95620,  0.45316),
        (97065,  0.46240), (99043,  0.47079), (100641, 0.47751),
        (101783, 0.48255), (102848, 0.48759), (103989, 0.49221),
        (106120, 0.49767),
    ])
 
    # PR vs corrected mass flow (Fig.19 – operating line okuma)
    pr_mf_raw = np.array([
        (0.28024, 1.6703), (0.31073, 1.8370), (0.33970, 2.0072),
        (0.36460, 2.1739), (0.38900, 2.3370), (0.40883, 2.4891),
        (0.42154, 2.5942), (0.43781, 2.7283), (0.45053, 2.8514),
        (0.45969, 2.9565), (0.46986, 3.0652), (0.47953, 3.1667),
        (0.48767, 3.2826), (0.49480, 3.3659), (0.50141, 3.4565),
        (0.50599, 3.5217), (0.51109, 3.6051), (0.51822, 3.7138),
    ])
 
    # η_c vs corrected mass flow (Fig.19 – verim kontur okuması)
    eta_mf_raw = np.array([
        (0.16117, 0.68875), (0.18204, 0.70483), (0.22330, 0.74505),
        (0.28058, 0.76890), (0.31019, 0.77805), (0.33932, 0.77666),
        (0.36505, 0.77444), (0.38835, 0.77583), (0.40825, 0.77278),
        (0.42087, 0.76613), (0.43689, 0.76474), (0.45049, 0.76002),
        (0.45922, 0.75697), (0.46942, 0.75475), (0.47864, 0.74754),
        (0.48738, 0.74754), (0.49417, 0.73922), (0.50097, 0.73368),
        (0.50583, 0.72813), (0.51068, 0.72397), (0.51748, 0.71787),
    ])
 
    # ── Interpolasyon ────────────────────────────────────────────────────────
    f_mfr  = interp1d(mfr_raw[:,0], mfr_raw[:,1], kind='cubic',
                      fill_value='extrapolate')
    f_pr   = interp1d(pr_mf_raw[:,0],  pr_mf_raw[:,1],  kind='cubic',
                      fill_value='extrapolate')
    f_etac = interp1d(eta_mf_raw[:,0], eta_mf_raw[:,1], kind='cubic',
                      fill_value='extrapolate')
 
    # ── Operating Line: N_corr → Nr → PR, η_c ───────────────────────────────
    N_grid   = np.linspace(62000, 106000, 300)
    Nr_grid  = N_grid / N_design_corr
    mf_grid  = f_mfr(N_grid)
    PR_grid  = f_pr(mf_grid)
    eta_grid = f_etac(mf_grid)

    # YENİ EKLENEN KISIM: Veriden mutlak zirveyi bulma (Göz kararı değil!)
    peak_idx   = np.argmax(eta_grid)      # Verimin en yüksek olduğu indeks
    eta_c_peak = eta_grid[peak_idx]       # Zirve verimi (0.778 civarı çıkacak)
    Nr_peak    = Nr_grid[peak_idx]        # Bu verime denk gelen devir oranı
 
    # ── Curve Fit ────────────────────────────────────────────────────────────
    def pr_model(Nr, a):
        return 1.0 + (PR_design - 1.0) * Nr**a
 
    def eta_model(Nr, kc):
        return eta_c_peak * np.exp(-kc * (Nr - Nr_peak)**2)
 
    popt_a,  _ = curve_fit(pr_model,  Nr_grid, PR_grid,  p0=[2.0],
                            bounds=(0.5, 5.0))
    popt_kc, _ = curve_fit(eta_model, Nr_grid, eta_grid, p0=[0.5],
                            bounds=(1e-4, 10.0))
 
    a_fit  = float(popt_a[0])
    kc_fit = float(popt_kc[0])
 
    # R² hesabı
    def r2(y_obs, y_pred):
        ss_res = np.sum((y_obs - y_pred)**2)
        ss_tot = np.sum((y_obs - y_obs.mean())**2)
        return 1.0 - ss_res / ss_tot
 
    R2_PR  = r2(PR_grid,  pr_model(Nr_grid,  a_fit))
    R2_eta = r2(eta_grid, eta_model(Nr_grid, kc_fit))
 
    return {
        "a_fit":    a_fit,
        "kc_fit":   kc_fit,
        "R2_PR":    R2_PR,
        "R2_eta":   R2_eta,
        "Nr_grid":  Nr_grid,
        "PR_grid":  PR_grid,
        "eta_grid": eta_grid,
        "eta_c_peak": float(eta_c_peak), # Numpy array'den float'a çeviriyoruz
        "Nr_peak":    float(Nr_peak)
    }

def synthetic_map_hybrid(
    Nr_corr, a, kc, kt,
    eta_c_peak, Nr_peak,
    PR_design=3.3812, eta_c_design=0.744, eta_t_design=0.850,
):
    PR    = max(1.01, 1.0 + (PR_design    - 1.0) * Nr_corr ** a)
    eta_c = float(np.clip(eta_c_peak * np.exp(-kc * (Nr_corr - Nr_peak)**2), 0.50, 0.95))
    eta_t = float(np.clip(eta_t_design * np.exp(-kt * (Nr_corr - 1.0)**2), 0.50, 0.95))
    return PR, eta_c, eta_t

def run_turbojet_with_synthetic_maps(
    h, M_ambient, rp_comp_des, T4_max, fuel_blend,
    N_ratio_mech,
    eta_c_peak, Nr_peak,
    eta_c_des,   
    eta_t_des,   
    a_exp,        
    kc,          
    kt,           
    eta_nozz=0.97, eta_mech=0.99, eta_comb=0.98, deltaP_comb=0.05
):
    # 1. Cantera Çözücü Hazırlığı
    gas = ct.Solution("air.yaml")
    T_sl, P_sl = 288.15, 101325.0
    g = 9.80665
    R = ct.gas_constant / gas.mean_molecular_weight

    # 2. ISA Atmosfer Modeli
    if h <= 11000:
        L = -0.0065
        T_static = T_sl + L * h
        P_static = P_sl * (T_static / T_sl) ** (-g / (L * R))
    else:
        T11 = T_sl - 0.0065 * 11000
        P11 = P_sl * (T11 / T_sl) ** (-g / (-0.0065 * R))
        T_static = T11
        P_static = P11 * np.exp(-g * (h - 11000) / (R * T11))

    # 3. Cantera ile Toplam (Stagnation) Durum Hesabı
    gas.TP = T_static, P_static
    h_static = gas.h
    V0 = M_ambient * gas.sound_speed
    
    # Enerji Dengesi: Toplam Entalpi = Statik Entalpi + Kinetik Enerji
    h_total = h_static + 0.5 * V0**2
    
    # Toplam entalpiden Toplam Sıcaklığı (T_tot) çekiyoruz
    # Basınç burada T_tot sonucunu etkilemez (Thermally Perfect Gas)
    gas.HP = h_total, P_static 
    T_tot = gas.T


    # --------------------- #ÖNEMLİ# -------------------------------
    Nr_corr = N_ratio_mech * np.sqrt(T_sl / T_tot)
 
    PR_comp, eta_comp_real, eta_turb_real = synthetic_map_hybrid(
        Nr_corr,
        PR_design=rp_comp_des,
        eta_c_design=eta_c_des,
        eta_t_design=eta_t_des,
        a=a_exp, kc=kc, kt=kt, 
        eta_c_peak=eta_c_peak, Nr_peak=Nr_peak
    )
 
    r = actual_brayton_cycle_turbojet(
        h=h, M_ambient=M_ambient, rp_comp=PR_comp,
        eta_comp=eta_comp_real, eta_turb=eta_turb_real, eta_nozz=eta_nozz,
        T4=T4_max, fuel_blend=fuel_blend,
        deltaP_comb=deltaP_comb, eta_mech=eta_mech, eta_comb=eta_comb
    )
 
    r["N_ratio_mech"] = N_ratio_mech
    r["N_corr"]       = Nr_corr
    r["PR_comp_used"] = PR_comp
    r["T4_used"]      = T4_max
    return r

if __name__ == "__main__":

    print("\n" + "★" * 72)
    print("  SAF-ThermoTwin |  Hibrit Motor Mimarisi: AMT Olympus HP")
    print("★" * 72)

    # =========================================================================
    # ADIM 1: AMT OLYMPUS HP TASARIM NOKTASI VERİLERİ
    # =========================================================================
    # Kaynak: Leylek (2012) Şekil 10, 18, 19 okuma değerleri
    # ─────────────────────────────────────────────────────────────────────────
    h_SLS         = 0.0
    M_SLS         = 0.0
    PR_design     = 3.3812   # Fig.19 operating line, 106,064 RPM
    EGT_leylek_K  = 977.63 # Fig.18 gerçek test
    eta_comp_dp   = 0.744  # Fig.19 verim haritası
    eta_turb_dp   = 0.85  # Leylek Fig.10, rp_turb≈1.72, N=59047 eğrisi
    eta_nozz_dp   = 0.97
    eta_mech_dp   = 0.99
    eta_comb_dp   = 0.98
    dP_comb_dp    = 0.05

    print(f"\n{'═'*72}")
    print("  ADIM 1: TASARIM NOKTASI VERİLERİ  (AMT Olympus HP @ SLS)")
    print(f"{'═'*72}")
    print(f"  PR_komp     = {PR_design:.2f}   (Leylek Tablo 1)")
    print(f"  EGT hedef   = {EGT_leylek_K:.1f} K  (Leylek Fig. 18, test verisi)")
    print(f"  η_komp_is   = {eta_comp_dp:.3f}  (Leylek Fig. 19)")
    print(f"  η_türbin_tt = {eta_turb_dp:.3f}  (Leylek Fig. 10 CFD)")

    # =========================================================================
    # ADIM 2: T4 KALİBRASYONU — EGT DOĞRULAMASI
    # =========================================================================
    # Amaç: Nozül çıkışı statik sıcaklığı (state6["T"]) ≈ EGT_leylek_K olacak
    # şekilde Türbin Giriş Sıcaklığı (T4) bulunur.
    # Yöntem: Biseksiyon (İkiye Bölme) algoritması
    # ─────────────────────────────────────────────────────────────────────────
    print(f"\n{'═'*72}")
    print("  ADIM 2: T4 KALİBRASYONU - EGT - JP4")
    print(f"{'═'*72}")
    print(f"  Hedef: T_nozül_çıkış (T6) = {EGT_leylek_K:.1f} K @ SLS, PR={PR_design}")
    print("  Yöntem: Biseksiyon – 60 iterasyon maks.\n")

    T4_lo, T4_hi = 900.0, 1800.0
    T4_cal = 0.5 * (T4_lo + T4_hi)
    converged = False

    for i_cal in range(60):
        T4_mid = 0.5 * (T4_lo + T4_hi)
        r_test = actual_brayton_cycle_turbojet(
            h_SLS, M_SLS, PR_design,
            eta_comp_dp, eta_turb_dp, eta_nozz_dp, T4_mid,
            compute_blend_properties_jp4("JP-4", 0.0, verbose=False),
            deltaP_comb=dP_comb_dp, eta_mech=eta_mech_dp, eta_comb=eta_comb_dp
        )
        T6_try = r_test["state6"]["T"]
        err    = T6_try - EGT_leylek_K
        if i_cal % 10 == 0:
            print(f"    iter={i_cal:3d}  T4={T4_mid:.2f} K  T6={T6_try:.2f} K  hata={err:+.3f} K")
        if abs(err) < 0.20:
            T4_cal    = T4_mid
            converged = True
            break
        if T6_try < EGT_leylek_K:
            T4_lo = T4_mid
        else:
            T4_hi = T4_mid
    else:
        T4_cal = T4_mid

    # Kalibre T4 ile son SLS hesabı
    r_sls = actual_brayton_cycle_turbojet(
        h_SLS, M_SLS, PR_design,
        eta_comp_dp, eta_turb_dp, eta_nozz_dp, T4_cal,
            compute_blend_properties_jp4("JP-4", 0.0, verbose=False),
        deltaP_comb=dP_comb_dp, eta_mech=eta_mech_dp, eta_comb=eta_comb_dp
    )
 
    EGT_pred  = r_sls["state6"]["T"]
    EGT_error = (EGT_pred - EGT_leylek_K) / EGT_leylek_K * 100
 
    print(f"\n  {'─'*50}")
    print(f"  T4_kalibre  = {T4_cal:.4f} K  ({'YAKINSADI ✔' if converged else 'yakınsamadı'})")
    print(f"  T6_tahmin   = {EGT_pred:.4f} K")
    print(f"  T6_hedef    = {EGT_leylek_K:.2f} K  (Leylek Fig.18)")
    print(f"  EGT hatası  = {EGT_error:+.4f} %")
    print(f"  FAR (SLS)   = {r_sls['FAR']:.5f}")
    print(f"  BWR (SLS)   = {r_sls['bwr']:.5f}")
    print(f"\n  {'─'*55}")
    print(f"  TABLE III — AMT OLYMPUS HP ENGINE STATE POINTS (SLS, JP-4)")
    print(f"  {'─'*55}")
    print(f"  {'Station':<30} {'T (K)':>10} {'P (kPa)':>10}")
    print(f"  {'─'*55}")
    print(f"  {'Ambient (SLS)':<30} {r_sls['state1']['T']:>10.2f} {r_sls['state1']['P']:>10.2f}")
    print(f"  {'Diffuser outlet':<30} {r_sls['state2']['T']:>10.2f} {r_sls['state2']['P']:>10.2f}")
    print(f"  {'Compressor outlet':<30} {r_sls['state3']['T']:>10.2f} {r_sls['state3']['P']:>10.2f}")
    print(f"  {'Turbine Inlet (T4)':<30} {r_sls['state4']['T']:>10.2f} {r_sls['state4']['P']:>10.2f}")
    print(f"  {'Turbine Outlet':<30} {r_sls['state5']['T']:>10.2f} {r_sls['state5']['P']:>10.2f}")
    print(f"  {'Nozzle Exit (EGT)':<30} {r_sls['state6']['T']:>10.2f} {r_sls['state6']['P']:>10.2f}")
    print(f"  {'─'*55}")
 
    # =========================================================================
    # ADIM 3 — SENTETİK HARİTA KALİBRASYONU (Leylek operating line)
    # =========================================================================
    print("\n" + "═" * 72)
    print("  ADIM 3: SENTETİK HARİTA KALİBRASYONU")
    print("  (Leylek Fig.16+19 sayısallaştırma verisi → curve_fit)")
    print("═" * 72)
 
    # calibrate_hybrid_map_from_leylek() → a, kc, R² değerlerini döndürür
    map_cal = calibrate_hybrid_map_from_leylek(
        PR_design=PR_design,
        eta_c_design=eta_comp_dp,
        
    )
    a_fit  = map_cal["a_fit"]
    kc_fit = map_cal["kc_fit"]
    eta_c_peak = map_cal["eta_c_peak"]
    Nr_peak    = map_cal["Nr_peak"]
 
    print(f"  Power-law üssü  a   = {a_fit:.4f}  (R²={map_cal['R2_PR']:.4f})")
    print(f"  Gaussian katsayı kc = {kc_fit:.4f}  (R²={map_cal['R2_eta']:.4f})")
    print(f"  Türbin sabiti   kt  = 0.20   (ampirik, radyal türbin)")
 
    # =========================================================================
    # ADIM 4 — SLS DOĞRULAMA (harita ile)
    # =========================================================================
    print("\n" + "═" * 72)
    print("  ADIM 4: SLS DOĞRULAMA")
    print("═" * 72)
 

    r_val = run_turbojet_with_synthetic_maps(
        h=h_SLS, M_ambient=M_SLS,
        rp_comp_des=PR_design, T4_max=T4_cal, fuel_blend=compute_blend_properties_jp4("JP-4", 0.0, verbose=False),
        N_ratio_mech=1.0, eta_c_peak=eta_c_peak, Nr_peak=Nr_peak,
        eta_c_des=eta_comp_dp, eta_t_des=eta_turb_dp,
        a_exp=a_fit, kc=kc_fit, kt=0.20,
        eta_nozz=eta_nozz_dp, eta_mech=eta_mech_dp,
        eta_comb=eta_comb_dp, deltaP_comb=dP_comb_dp
    )
 
    EGT_map_err = (r_val["state6"]["T"] - EGT_leylek_K) / EGT_leylek_K * 100
    print(f"  Nr_corr (SLS) = {r_val['N_corr']:.4f}  (beklenen: 1.0000)")
    print(f"  PR_comp       = {r_val['PR_comp_used']:.4f}  (beklenen: {PR_design})")
    print(f"  η_c           = {r_val['eta_comp_used']:.4f}")
    print(f"  η_t           = {r_val['eta_turb_used']:.4f}")
    print(f"  T6/EGT        = {r_val['state6']['T']:.4f} K")
    print(f"  EGT hatası    = {EGT_map_err:+.4f} %  (Leylek: {EGT_leylek_K} K)")

    # =========================================================================
    # ADIM 5.a — CRUISE BASELINE (10,000 m, Mach 0.6, JP-4)
    # =========================================================================
    print("\n" + "═" * 72)
    print("  ADIM 5.a: CRUISE BASELINE (10,000 m, M=0.6, JP-4)")
    print("═" * 72)
 
    # N_ratio_mech=1.0 → motor aynı fiziksel devirde çalışıyor
    # Daha soğuk giriş (T_tot=239.4 K) → Nr_corr = sqrt(288.15/239.40) = 1.097
    r_cruise_jp4 = run_turbojet_with_synthetic_maps(
        h=10000, M_ambient=0.6,
        rp_comp_des=PR_design, T4_max=T4_cal, fuel_blend=compute_blend_properties_jp4("JP-4", 0.0, verbose=False),
        N_ratio_mech=1.0,
        eta_c_peak=eta_c_peak, Nr_peak=Nr_peak,
        eta_c_des=eta_comp_dp, eta_t_des=eta_turb_dp,
        a_exp=a_fit, kc=kc_fit, kt=0.20,
        eta_nozz=eta_nozz_dp, eta_mech=eta_mech_dp,
        eta_comb=eta_comb_dp, deltaP_comb=dP_comb_dp
    )

    print(f"  Nr_corr       = {r_cruise_jp4['N_corr']:.4f}  (aynı fiziksel RPM, düşük T_tot)")
    print(f"  PR_comp       = {r_cruise_jp4['PR_comp_used']:.4f}")
    print(f"  η_c           = {r_cruise_jp4['eta_comp_used']:.4f}")
    print(f"  η_t           = {r_cruise_jp4['eta_turb_used']:.4f}")
    print(f"  T3 (K)        = {r_cruise_jp4['state3']['T']:.2f}")
    print(f"  T4 (K)        = {r_cruise_jp4['T4_used']:.2f}  (sabit)")
    print(f"  T5 (K)        = {r_cruise_jp4['state5']['T']:.2f}")
    print(f"  T6/EGT (K)    = {r_cruise_jp4['state6']['T']:.2f}")
    print(f"  Fs [N/(kg/s)] = {r_cruise_jp4['specific_thrust']:.2f}")
    print(f"  SFC [mg/(N·s)]= {r_cruise_jp4['SFC_mg_Ns']:.4f}")
    print(f"  FAR           = {r_cruise_jp4['FAR']:.6f}")
    print(f"  η_th (%)      = {r_cruise_jp4['thermal_efficiency']*100:.2f}")
    print(f"  η_p  (%)      = {r_cruise_jp4['propulsive_efficiency']*100:.2f}")
    print(f"  η_o  (%)      = {r_cruise_jp4['overall_efficiency']*100:.2f}")
    print(f"  CO2_net (g/kg)= {r_cruise_jp4['net_CO2_per_kgair']*1000:.2f}")
    print(f"  Choked?       = {r_cruise_jp4['choked']}")

    # =========================================================================
    # ADIM 5.b — CRUISE BASELINE (10,000 m, Mach 0.6, Jet-A)
    # =========================================================================
    print("\n" + "═" * 72)
    print("  ADIM 5.b: CRUISE BASELINE (10,000 m, M=0.6, Jet-A)")
    print("═" * 72)
 
    # N_ratio_mech=1.0 → motor aynı fiziksel devirde çalışıyor
    # Daha soğuk giriş (T_tot=239.4 K) → Nr_corr = sqrt(288.15/239.40) = 1.097
    r_cruise = run_turbojet_with_synthetic_maps(
        h=10000, M_ambient=0.6,
        rp_comp_des=PR_design, T4_max=T4_cal, fuel_blend=compute_blend_properties("Jet-A", 0.0, verbose=False),
        N_ratio_mech=1.0,
        eta_c_peak=eta_c_peak, Nr_peak=Nr_peak,
        eta_c_des=eta_comp_dp, eta_t_des=eta_turb_dp,
        a_exp=a_fit, kc=kc_fit, kt=0.20,
        eta_nozz=eta_nozz_dp, eta_mech=eta_mech_dp,
        eta_comb=eta_comb_dp, deltaP_comb=dP_comb_dp
    )

    print(f"  Nr_corr       = {r_cruise['N_corr']:.4f}  (aynı fiziksel RPM, düşük T_tot)")
    print(f"  PR_comp       = {r_cruise['PR_comp_used']:.4f}")
    print(f"  η_c           = {r_cruise['eta_comp_used']:.4f}")
    print(f"  η_t           = {r_cruise['eta_turb_used']:.4f}")
    print(f"  T3 (K)        = {r_cruise['state3']['T']:.2f}")
    print(f"  T4 (K)        = {r_cruise['T4_used']:.2f}  (sabit)")
    print(f"  T5 (K)        = {r_cruise['state5']['T']:.2f}")
    print(f"  T6/EGT (K)    = {r_cruise['state6']['T']:.2f}")
    print(f"  Fs [N/(kg/s)] = {r_cruise['specific_thrust']:.2f}")
    print(f"  SFC [mg/(N·s)]= {r_cruise['SFC_mg_Ns']:.4f}")
    print(f"  FAR           = {r_cruise['FAR']:.6f}")
    print(f"  η_th (%)      = {r_cruise['thermal_efficiency']*100:.2f}")
    print(f"  η_p  (%)      = {r_cruise['propulsive_efficiency']*100:.2f}")
    print(f"  η_o  (%)      = {r_cruise['overall_efficiency']*100:.2f}")
    print(f"  CO2_net (g/kg)= {r_cruise['net_CO2_per_kgair']*1000:.2f}")
    print(f"  Choked?       = {r_cruise['choked']}")

    # =========================================================================
    # ADIM 6 — SAF KARŞILAŞTIRMASI (10,000 m, Mach 0.6)
    # Maks. ASTM D7566 karışım fraksiyonlarında
    # =========================================================================

    print("\n" + "═" * 72)
    print("  ADIM 6: SAF KARŞILAŞTIRMASI (10,000 m, M=0.6)")
    print("  Maks. ASTM D7566 karışım fraksiyonlarında")
    print("═" * 72)
 
    # (saf_adı, karışım_fraksiyonu) — ASTM D7566 maks. sınırlar
    SAF_CASES = [
        ("HEFA-SPK",      0.50),
        ("FT-SPK",        0.50),
        ("ATJ-SPK",       0.50),
        ("HFS-SIP",       0.10),
        ("Co-processing", 0.50),
    ]

    # Baseline referans değerleri (Jet-A, Adım 5'ten)
    Fs_base   = r_cruise["specific_thrust"]
    SFC_base  = r_cruise["SFC_mg_Ns"]
    LCA_base  = compute_blend_properties("Jet-A", 0.0, verbose=False)["lifecycle_CO2_gMJ"]  # 89.0

    # Tablo başlığı
    print(f"\n  {'Yakıt':<22} {'Blend':>5}  {'Fs':>9}  {'ΔSFC(%)':>9}  "
          f"{'EGT(K)':>8}  {'CO2_net(g/kg)':>14}  {'LCA tasarrufu(%)':>17}")
    print(f"  {'─'*22}  {'─'*5}  {'─'*9}  {'─'*9}  {'─'*8}  {'─'*14}  {'─'*17}")
 
    # Baseline satırı
    print(f"  {'Jet-A (baseline)':<22}  {'%0':>5}  "
          f"{Fs_base:>9.2f}  {'—':>9}  "
          f"{r_cruise['state6']['T']:>8.2f}  "
          f"{r_cruise['net_CO2_per_kgair']*1000:>14.2f}  "
          f"{'—':>17}")
    
    # Her SAF için döngü
    for saf_adi, frac in SAF_CASES:
        saf_blend = compute_blend_properties(saf_adi, frac, verbose=False)
 
        r_saf = run_turbojet_with_synthetic_maps(
            h=10000, M_ambient=0.6,
            rp_comp_des=PR_design, T4_max=T4_cal, fuel_blend=saf_blend,
            N_ratio_mech=1.0,
            eta_c_peak=eta_c_peak, Nr_peak=Nr_peak,
            eta_c_des=eta_comp_dp, eta_t_des=eta_turb_dp,
            a_exp=a_fit, kc=kc_fit, kt=0.20,
            eta_nozz=eta_nozz_dp, eta_mech=eta_mech_dp,
            eta_comb=eta_comb_dp, deltaP_comb=dP_comb_dp
        )
 
        dSFC    = (r_saf["SFC_mg_Ns"] - SFC_base) / SFC_base * 100
        LCA_sav = (1.0 - saf_blend["lifecycle_CO2_gMJ"] / LCA_base) * 100
        etiket  = f"{saf_adi} (%{frac*100:.0f})"
 
        print(f"  {etiket:<22}  {frac*100:>4.0f}%  "
              f"{r_saf['specific_thrust']:>9.2f}  "
              f"{dSFC:>+9.3f}  "
              f"{r_saf['state6']['T']:>8.2f}  "
              f"{r_saf['net_CO2_per_kgair']*1000:>14.2f}  "
              f"{LCA_sav:>17.1f}")
 
    print(f"\n  {'═'*72}")
    print("  Analiz tamamlandı.")
    print(f"  {'═'*72}")