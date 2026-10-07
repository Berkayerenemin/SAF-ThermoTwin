# SAF-ThermoTwin

**Validated non-ideal Brayton cycle framework for thermodynamic performance and CORSIA-aligned lifecycle CO₂ analysis of Sustainable Aviation Fuel (SAF) blends in micro-turbojet engines.**

SAF-ThermoTwin is an open-source Python tool that couples a thermally perfect gas turbojet cycle model (resolved with [Cantera](https://cantera.org/)) with CORSIA-aligned lifecycle CO₂ accounting. It is calibrated against published sea-level static test data of the **AMT Olympus HP** micro-turbojet and used to compare ASTM D7566-certified SAF pathways at representative cruise conditions.

> **Paper:**
> B. E. Emin, *"SAF-ThermoTwin: Validated Non-Ideal Brayton Cycle Framework for Thermodynamic Performance and Emissions Modelling of SAF Blends in Micro Turbojet Engines,"* 2026 11th Int. Conf. on Recent Advances in Air and Space Technologies (RAST), Istanbul, Turkiye, 2026, pp. 1–6.
> DOI: [10.1109/RAST69551.2026.11672525](https://doi.org/10.1109/RAST69551.2026.11672525)
>
> **Author:** Berkay Eren Emin — Dept. of Mechanical Engineering, İzmir Katip Çelebi University
> ORCID: [0009-0000-1779-9443](https://orcid.org/0009-0000-1779-9443)

> [!NOTE]
> Console output and inline comments are currently in **Turkish**. A Turkish version of this README follows below.

---

## Key Features

- **Non-ideal Brayton cycle** with five components: inlet/diffuser, compressor, combustor, turbine, convergent nozzle.
- **Thermally perfect gas** thermodynamics (NASA polynomials via Cantera): temperature-dependent cₚ, enthalpy and entropy; no constant-γ shortcuts (nozzle exit velocity is computed from the enthalpy drop).
- **Synthetic off-design performance map** fitted to the digitised compressor operating line from Leylek (2012):
  - Power-law pressure ratio: `PRc = 1 + (PRc,design − 1) · Nr^a`
  - Gaussian efficiency decay for compressor and turbine
- **Design-point calibration:** turbine inlet temperature (T4) is back-calculated by bisection from the measured EGT (choked nozzle at SLS, JP-4).
- **Fuel database and blending model:** Jet-A, HEFA-SPK, FT-SPK, ATJ-SPK, HFS-SIP, with volume→mass blend conversion, mass-weighted LHV, ideal-mixing density and ASTM D7566 maximum blend limits.
- **CORSIA-aligned lifecycle CO₂** savings and biogenic-corrected net CO₂ per kg of air.
- **Fixed-control strategy** (constant shaft speed and T4): differences between fuels arise only from fuel properties.

## Model Overview

| Step | Description |
|------|-------------|
| 1 | Design-point data (AMT Olympus HP @ SLS) |
| 2 | T4 calibration via bisection against measured EGT (JP-4, SLS) |
| 3 | Synthetic performance map calibration (nonlinear regression on digitised Leylek map data) |
| 4 | SLS validation with synthetic maps |
| 5a | Cruise baseline (10,000 m, Mach 0.6, JP-4) |
| 5b | Cruise baseline (10,000 m, Mach 0.6, Jet-A) |
| 6 | Multi-SAF comparison at maximum ASTM D7566 blend fractions |

Fuel-to-air ratio at fixed T4:

```
f = (h4 − h3) / (ηb · LHV_blend)
```

Specific thrust (with fuel mass addition):

```
Fs = (1 + f)·V6 − V0 + (P6 − P1)(1 + f) / (ρ6·V6)
```

## Baseline Engine (AMT Olympus HP, JP-4, SLS)

| Parameter | Value |
|-----------|-------|
| Compressor pressure ratio | 3.3812 |
| Design-point compressor η_is | 0.744 |
| Design-point turbine η_is | 0.85 |
| Calibrated T4 | 1264.31 K |
| Combustor efficiency / pressure drop | 0.98 / 5 % |
| Mechanical efficiency | 0.99 |
| Nozzle efficiency | 0.97 |
| Design speed | 106,064 RPM |

## Validation

Reference: Z. Leylek, *An investigation into performance modelling of a small gas turbine engine*, DSTO-TR-2757, 2012.

| Metric | Reference | Model | Error |
|--------|-----------|-------|-------|
| EGT (SLS, 106,064 RPM) | 977.63 K | 977.66 K | +0.003 % |

Map fit quality: `a = 2.4089` (R² = 0.9960), `kc = 0.4871` (R² = 0.9325); turbine decay `kt = 0.20` (fixed).

> [!NOTE]
> T4 is not directly measured in the reference data; it is calibrated so that the model reproduces the measured EGT. The off-design maps are validated through physical consistency with the published compressor map data.

## Results (10,000 m, Mach 0.6, Maximum ASTM D7566 Blend)

| Fuel (blend) | Fs [N/(kg/s)] | ΔSFC [%] | EGT [K] | CO₂,net [g/kg air] | Lifecycle CO₂ saving [%] |
|--------------|:-------------:|:--------:|:-------:|:-------------------:|:------------------------:|
| Jet-A (0 %)     | 580.87 | —      | 966.28 | 70.69 | —    |
| HEFA-SPK (50 %) | 580.65 | −0.966 | 966.26 | 35.73 | 36.4 |
| FT-SPK (50 %)   | 580.60 | −1.163 | 966.25 | 35.96 | 42.7 |
| ATJ-SPK (50 %)  | 580.70 | −0.751 | 966.26 | 35.67 | 30.6 |
| HFS-SIP (10 %)  | 580.85 | −0.086 | 966.28 | 63.75 | 4.9  |

Jet-A cruise baseline: SFC = 38.65 mg/(N·s), η_th = 17.76 %, η_p = 60.69 %, η_o = 10.78 %.

## Installation

Requirements: Python 3.10+ and the packages in `requirements.txt`.

```bash
git clone https://github.com/<your-username>/SAF-ThermoTwin.git
cd SAF-ThermoTwin
python -m venv .venv
source .venv/bin/activate        # Windows: .venv\Scripts\activate
pip install -r requirements.txt
```

## Usage

Run the full workflow (calibration → map fitting → SLS validation → cruise baseline → SAF comparison):

```bash
python saf_thermotwin.py
```

To explore other cases, edit the cruise condition (`h`, `M_ambient`) or the `SAF_CASES` list near the end of the script, and extend `FUEL_DATABASE` with new fuel entries (H/C ratio, LHV, density, CO₂ factor, biogenic carbon fraction, lifecycle CO₂ intensity).

## Repository Structure

```
.
├── saf_thermotwin.py      # Full framework (fuel DB, cycle model, maps, workflow)
├── requirements.txt
├── CITATION.cff
├── LICENSE                # GNU GPL v3.0
└── README.md              # This file (English + Turkish)
```

> **Note on filename:** The script was previously named `SAF-ThermoTwin-AMTOlympusHP-SON-Rev2.py`. It has been renamed to `saf_thermotwin.py` for importability and readability; the content is identical.

## Assumptions and Limitations

- Single-spool, afterburner-free micro-turbojet (AMT Olympus HP); other architectures require re-calibration.
- Engine runs at fixed maximum continuous mechanical speed and fixed calibrated T4; fuel-to-air ratio adjusts passively.
- Design-point efficiencies come from digitised published map data; digitisation uncertainty is not propagated.
- Lifecycle CO₂ values are literature medians (CORSIA-based) and vary widely with feedstock and region (e.g. HEFA: roughly 10–88 gCO₂e/MJ). Results are indicative, not certified LCA values.
- Fuels are represented by surrogate compositions with mass-weighted linear blending rules.
- The code contains a **"Co-processing"** entry that is **not** an ASTM D7566 pathway and is **not** part of the paper's results; its lifecycle value is indicative only.
- No combustion chemistry, non-CO₂ effects (NOₓ, soot, contrails) or transient behaviour is modelled.

## Data Sources

The synthetic performance map is built from digitised operating-line data from Leylek (2012) DSTO-TR-2757. The original report PDF and figures are not included in this repository; only the derived numerical arrays embedded in `calibrate_hybrid_map_from_leylek()` are used.

## Citation

If you use this software, please cite both the code (see `CITATION.cff`) and the associated paper:

```bibtex
@inproceedings{emin2026saftt,
  author    = {Emin, Berkay Eren},
  title     = {SAF-ThermoTwin: Validated Non-Ideal Brayton Cycle Framework
               for Thermodynamic Performance and Emissions Modelling of SAF
               Blends in Micro Turbojet Engines},
  booktitle = {2026 11th International Conference on Recent Advances in Air
               and Space Technologies (RAST)},
  year      = {2026},
  pages     = {1--6},
  address   = {Istanbul, Turkiye},
  publisher = {IEEE},
  doi       = {10.1109/RAST69551.2026.11672525}
}
```

## References

1. Z. Leylek, *An investigation into performance modelling of a small gas turbine engine*, DSTO-TR-2757, 2012.
2. M. A. Aksoy and O. Son, "Performance modelling and efficiency analysis of a micro turbojet engine under varying flight conditions", *J. Braz. Soc. Mech. Sci. Eng.*, 48, 238, 2026.
3. M. Prussi et al., "CORSIA: The first internationally adopted approach to calculate lifecycle GHG emissions for aviation fuels", *Renew. Sustain. Energy Rev.*, 150, 111398, 2021.
4. ASTM D7566-23, *Standard Specification for Aviation Turbine Fuel Containing Synthesized Hydrocarbons*.
5. D. S. Lee et al., "The contribution of global aviation to anthropogenic climate forcing for 2000 to 2018", *Atmos. Environ.*, 244, 117834, 2021.

## License

Released under the **GNU General Public License v3.0** — see [LICENSE](LICENSE) for details.

GPL v3.0 was chosen to require that derivative works remain equally transparent and openly accessible to the academic community, enabling independent verification and extension of the framework.

## Disclaimer

This is a research tool. Outputs are model predictions and must not be used for flight-critical, certification or regulatory purposes without independent validation.

---
---

# SAF-ThermoTwin — Türkçe / Turkish

**AMT Olympus HP mikro-turbojet motorunda SAF karışımlarının termodinamik performansı ve CORSIA uyumlu yaşam döngüsü CO₂ analizinin doğrulanmış ideal-olmayan Brayton çevrimi çerçevesi.**

SAF-ThermoTwin, [Cantera](https://cantera.org/) tabanlı ısıl-mükemmel gaz termodinamiğini CORSIA uyumlu yaşam döngüsü CO₂ muhasebesiyle birleştiren açık kaynaklı bir Python aracıdır. **AMT Olympus HP** mikro-turbojet motorunun deniz seviyesi statik test verilerine göre kalibre edilmiş; temsili seyir koşullarında ASTM D7566 onaylı SAF yollarını karşılaştırmak için kullanılmıştır.

> **Makale:**
> B. E. Emin, *"SAF-ThermoTwin: Validated Non-Ideal Brayton Cycle Framework for Thermodynamic Performance and Emissions Modelling of SAF Blends in Micro Turbojet Engines,"* 2026 11th Int. Conf. on Recent Advances in Air and Space Technologies (RAST), İstanbul, Türkiye, 2026, pp. 1–6.
> DOI: [10.1109/RAST69551.2026.11672525](https://doi.org/10.1109/RAST69551.2026.11672525)
>
> **Yazar:** Berkay Eren Emin — Makine Mühendisliği Bölümü, İzmir Katip Çelebi Üniversitesi
> ORCID: [0009-0000-1779-9443](https://orcid.org/0009-0000-1779-9443)

## Temel Özellikler

- **İdeal-olmayan Brayton çevrimi** — beş bileşen: hava alığı/difüzör, kompresör, yanma odası, türbin, yakınsak nozul.
- **Isıl-mükemmel gaz termodinamiği** — Cantera üzerinden NASA polinom: sıcaklık bağımlı cₚ, entalpi ve entropi; sabit γ yaklaşımı yok; nozul çıkış hızı entalpi düşüşünden hesaplanır.
- **Sentetik off-design performans haritası** — Leylek (2012)'den sayısallaştırılan işletme doğrusuna nonlineer regresyon:
  - Üs-kuvveti basınç oranı: `PRc = 1 + (PRc,tasarım − 1) · Nr^a`
  - Kompresör ve türbin verimi için Gaussian bozunumu
- **Tasarım noktası kalibrasyonu:** T4, ölçülen EGT'ye göre biseksiyon yöntemiyle geri hesaplanır (boğulmuş nozul, SLS, JP-4).
- **Yakıt veritabanı ve karışım modeli:** Jet-A, HEFA-SPK, FT-SPK, ATJ-SPK, HFS-SIP; hacim→kütle dönüşümü, kütleye göre LHV, ideal karışım yoğunluğu ve ASTM D7566 maks. karışım sınırları.
- **CORSIA uyumlu yaşam döngüsü CO₂** tasarrufu ve kg hava başına biyojenik-düzeltilmiş net CO₂.
- **Sabit kontrol stratejisi** (sabit şaft devri ve T4) — yakıtlar arası farklar yalnızca yakıt özelliklerinden kaynaklanır.

## Model Özeti

| Adım | Açıklama |
|------|----------|
| 1 | Tasarım-noktası verileri (AMT Olympus HP @ SLS) |
| 2 | T4 kalibrasyonu — ölçülen EGT'ye biseksiyon (JP-4, SLS) |
| 3 | Sentetik harita kalibrasyonu (Leylek işletme doğrusu, nonlineer regresyon) |
| 4 | SLS doğrulaması (sentetik haritayla) |
| 5a | Seyir bazlı (10.000 m, Mach 0,6, JP-4) |
| 5b | Seyir bazlı (10.000 m, Mach 0,6, Jet-A) |
| 6 | Maks. ASTM D7566 karışımlarında çoklu SAF karşılaştırması |

Sabit T4'te yakıt-hava oranı:

```
f = (h4 − h3) / (ηb · LHV_karışım)
```

Özgül itki (yakıt kütlesi katkısıyla):

```
Fs = (1 + f)·V6 − V0 + (P6 − P1)(1 + f) / (ρ6·V6)
```

## Referans Motor (AMT Olympus HP, JP-4, SLS)

| Parametre | Değer |
|-----------|-------|
| Kompresör basınç oranı | 3,3812 |
| Tasarım-noktası kompresör η_is | 0,744 |
| Tasarım-noktası türbin η_is | 0,85 |
| Kalibre T4 | 1264,31 K |
| Yanma odası verimi / basınç kaybı | 0,98 / %5 |
| Mekanik verim | 0,99 |
| Nozul verimi | 0,97 |
| Tasarım devri | 106.064 RPM |

## Doğrulama

Kaynak: Z. Leylek, *An investigation into performance modelling of a small gas turbine engine*, DSTO-TR-2757, 2012.

| Metrik | Referans | Model | Hata |
|--------|----------|-------|------|
| EGT (SLS, 106.064 RPM) | 977,63 K | 977,66 K | +%0,003 |

Harita uyum kalitesi: `a = 2,4089` (R² = 0,9960), `kc = 0,4871` (R² = 0,9325); türbin bozunumu `kt = 0,20` (sabit).

> [!NOTE]
> T4 referans veride doğrudan ölçülmez; model ölçülen EGT'yi yeniden üretecek şekilde kalibre edilir. Off-design haritaları, yayımlanan kompresör harita verisiyle fiziksel tutarlılık üzerinden doğrulanmıştır.

## Sonuçlar (10.000 m, Mach 0,6, Maks. ASTM D7566 Karışımı)

| Yakıt (karışım) | Fs [N/(kg/s)] | ΔSFC [%] | EGT [K] | CO₂,net [g/kg hava] | Yaşam döngüsü CO₂ tasarrufu [%] |
|-----------------|:-------------:|:--------:|:-------:|:--------------------:|:-------------------------------:|
| Jet-A (%0)       | 580,87 | —      | 966,28 | 70,69 | —    |
| HEFA-SPK (%50)  | 580,65 | −0,966 | 966,26 | 35,73 | 36,4 |
| FT-SPK (%50)    | 580,60 | −1,163 | 966,25 | 35,96 | 42,7 |
| ATJ-SPK (%50)   | 580,70 | −0,751 | 966,26 | 35,67 | 30,6 |
| HFS-SIP (%10)   | 580,85 | −0,086 | 966,28 | 63,75 | 4,9  |

Seyir bazlı Jet-A değerleri: SFC = 38,65 mg/(N·s), η_th = %17,76, η_p = %60,69, η_o = %10,78.

## Kurulum

Gereksinimler: Python 3.10+ ve `requirements.txt` içindeki paketler.

```bash
git clone https://github.com/<kullanıcı-adınız>/SAF-ThermoTwin.git
cd SAF-ThermoTwin
python -m venv .venv
.venv\Scripts\activate          # Linux/macOS: source .venv/bin/activate
pip install -r requirements.txt
```

## Kullanım

Tam iş akışını çalıştırın (kalibrasyon → harita uyumlama → SLS doğrulama → seyir bazlı → SAF karşılaştırması):

```bash
python saf_thermotwin.py
```

> [!NOTE]
> Konsol çıktıları Türkçe'dir.

Farklı koşullar için kodun sonundaki `h`, `M_ambient` değerlerini veya `SAF_CASES` listesini düzenleyin. Yeni yakıtlar için `FUEL_DATABASE`'i genişletin (H/C oranı, LHV, yoğunluk, CO₂ faktörü, biyojenik karbon fraksiyonu, yaşam döngüsü CO₂ yoğunluğu).

## Depo Yapısı

```
.
├── saf_thermotwin.py      # Tam çerçeve (yakıt VT, çevrim modeli, haritalar, iş akışı)
├── requirements.txt
├── CITATION.cff
├── LICENSE                # GNU GPL v3.0
└── README.md              # Bu dosya (İngilizce + Türkçe)
```

> **Dosya adı hakkında:** Orijinal adı `SAF-ThermoTwin-AMTOlympusHP-SON-Rev2.py` olan betik, import edilebilirlik ve okunabilirlik için `saf_thermotwin.py` olarak yeniden adlandırılmıştır; içerik aynıdır.

## Varsayımlar ve Kısıtlamalar

- Tek bobinli, ardyakıcısız mikro-turbojet (AMT Olympus HP); diğer mimariler yeniden kalibrasyon gerektirir.
- Motor sabit maks. mekanik devirde ve sabit kalibre T4'te çalışır; yakıt-hava oranı pasif olarak ayarlanır.
- Tasarım noktası verimleri sayısallaştırılmış yayın haritası verilerinden elde edilir; sayısallaştırma belirsizliği hesaplara yayılmaz.
- Yaşam döngüsü CO₂ değerleri literatür medyanlarıdır (CORSIA tabanlı) ve hammadde ile bölgeye göre önemli ölçüde değişir (ör. HEFA: yakl. 10–88 gCO₂e/MJ). Sonuçlar belirticidir, sertifikalı YDA değerleri değildir.
- Yakıtlar vekil bileşimlerle temsil edilir; kütleye göre doğrusal karışım kuralı uygulanır.
- **"Co-processing"** girişi ASTM D7566 kapsamında bir yol **değildir** ve **makale sonuçlarına dahil edilmemiştir**; yaşam döngüsü değeri yalnızca göstergeseldir.
- Yanma kimyası, CO₂-dışı etkiler (NOₓ, is, kondenstras) veya geçici davranış modellenmemiştir.

## Veri Kaynakları Hakkında

Sentetik performans haritası, **Leylek (2012) DSTO-TR-2757** raporundan sayısallaştırılmış grafik verisi kullanılarak oluşturulmuştur. Raporun PDF'i veya orijinal görseller repoya eklenmemiştir; yalnızca `calibrate_hybrid_map_from_leylek()` işlevi içindeki türetilmiş sayısal diziler kullanılmaktadır.

## Atıf

Bu yazılımı kullanıyorsanız lütfen hem kodu (`CITATION.cff`) hem de ilgili makaleyi atıf yapın (yukarıdaki BibTeX girişine bakın).

## Lisans

**GNU Genel Kamu Lisansı v3.0** kapsamında yayınlanmıştır — ayrıntılar için [LICENSE](LICENSE) dosyasına bakın.

GPL v3.0 seçilmesinin gerekçesi: Türev çalışmaların da aynı şeffaflıkla paylaşılması zorunlu kılınarak akademik topluluğun çalışmayı bağımsız olarak doğrulayabilmesi ve genişletebilmesi hedeflenmiştir.

## Sorumluluk Reddi

Bu bir araştırma aracıdır. Çıktılar model öngörüleridir; bağımsız doğrulama yapılmadan uçuş-kritik, sertifikasyon veya düzenleyici amaçlarla kullanılmamalıdır.
