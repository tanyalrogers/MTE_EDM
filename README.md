# Integration of Metabolic Theory of Ecology and Empirical Dynamic Modeling (MTE-EDM)

<!-- badges: start -->
<!-- badges: end -->

This repository contains the data and code associated with the paper "Constraining nonlinear time series modeling with the metabolic theory of ecology" by Stephan Munch, Tanya Rogers, Celia Symons, David Anderson, and Frank Pennekamp published in *PNAS* (2023).

![flowchart](flowchart.png)

Laboratory data are located in `data/constant_temp_expts/`, and field data files are located in `data/field_data_interpolated/`. For the field data, time is in units of decimal years, temperature in degrees C. See `data/MTE_EDM_metadata.csv` for the species information corresponding to each file. To comply with data licensing polices, the Lake Geneva data are omitted. The file `run_metadata_public.csv` can be substituted for `run_metadata.csv` to run the model-fitting code without these series. 

Data sources are listed below and in Table S2 of the publication. Reuse of the **source data** in this repository should cite both the original source(s) and this study, and comply with the original license terms for each dataset. The **code and derived outputs** in this repository (with the exception of the cycle period analysis) were developed by U.S. government employees (SM, TR) and are considered public domain. DA and TR coauthored the code for cycle period analysis; DA chooses to also make this code public domain.

### Lab data sources

DeLong, J.P.; Lyon, S. (2020). Temperature alters the shape of predator–prey cycles through effects on underlying mechanisms. PeerJ 8, e9377

Fussmann, K.E.; Schwarzmüller, F.; Brose, U.; Jousset, A.; Rall, B.C. (2014) Ecological stability in response to warming. Nat. Clim. Change 4, 206–210 

Halbach, U. (1984) Population dynamics of rotifers and its consequences for ecotoxicology. Hydrobiologia 109, 79–96

Laughton, A.M.; Knell, R.J. (2019) Warming at the population level: Effects on age structure, density, and generation cycles. Ecol. Evol. 9, 4403–4420 

### Field data sources

**Narragansett Bay**

Smayda, T.J. & the Bunker C community (1959-1997). Narragansett Bay Plankton Time Series. Graduate School of Oceanography, URI. Data available at: https://www.nabats.org/

**Lake Greifensee**

Thomas, M.K.; Fontana, S.; Reyes, M.; Kehoe, M.; Pomati, F. (2019). Data from: The predictability of a lake phytoplankton community, over time-scales of hours to years, Dryad, Dataset, https://doi.org/10.5061/dryad.r4454 

**Lake Geneva**

Rimet F, Anneville O, Barbet D, Chardon C, Crépin L, Domaizon I, Dorioz J-M, Espinat L, Frossard V, Guillard J, Goulon C, Hamelet V, Hustache J-C, Jacquet S, Lainé L, Montuelle B, Perney P, Quetin P, Rasconi S, Schellenberger A, Tran-Khac V, Monet G. (2020) The Observatory on LAkes (OLA) database: Sixty years of environmental data accessible to the public. J Limnol. 79. https://doi.org/10.4081/jlimnol.2020.1944 © SOERE OLA-IS, AnaEE-France, INRA of Thonon-les-Bains, CIPEL, Dec 19 2019, developed by the Eco-Informatique ORE system of the INRA.

**Tea Tortrix Moth - Japan**

Nelson, W.A.; Bjornstad, O.N.; Yamanaka, T. (2013). Data from: Recurrent insect outbreaks caused by temperature-driven changes in system stability, Dryad, Dataset, https://doi.org/10.5061/dryad.n11d4

**Wadden Sea**

Martens, P. (2007). Abundance of zooplankton at times series station List Reede in 1999-2006. Alfred Wegener Institute - Wadden Sea Station Sylt, PANGAEA, https://doi.org/10.1594/PANGAEA.646280; https://doi.org/10.1594/PANGAEA.646281; https://doi.org/10.1594/PANGAEA.646282; https://doi.org/10.1594/PANGAEA.646283; https://doi.org/10.1594/PANGAEA.646284; https://doi.org/10.1594/PANGAEA.646285; https://doi.org/10.1594/PANGAEA.646286; https://doi.org/10.1594/PANGAEA.646287

Martens, P. (2011). Abundance of zooplankton at times series station List Reede in 2007-2008. Alfred Wegener Institute - Wadden Sea Station Sylt, PANGAEA, https://doi.org/10.1594/PANGAEA.756061; https://doi.org/10.1594/PANGAEA.756062

van Beusekom, J. (2010). Hydrochemistry time series at station List Reede in 1999-2007. Alfred Wegener Institute - Wadden Sea Station Sylt, PANGAEA 
https://doi.org/10.1594/PANGAEA.745142; https://doi.org/10.1594/PANGAEA.745143; https://doi.org/10.1594/PANGAEA.745144; https://doi.org/10.1594/PANGAEA.745145; https://doi.org/10.1594/PANGAEA.745146; https://doi.org/10.1594/PANGAEA.745147; https://doi.org/10.1594/PANGAEA.745148; https://doi.org/10.1594/PANGAEA.745149; https://doi.org/10.1594/PANGAEA.745150; https://doi.org/10.1594/PANGAEA.745151

**Greece**

Damos, P. (2016). Using multivariate cross correlations, Granger causality and graphical models to quantify spatiotemporal synchronization and causality between pest populations. BMC Ecology, 16, 33. https://doi.org/10.1186/s12898-016-0087-7

**Bermuda**

Bermuda Atlantic Time Series: http://batsftp.bios.edu/BATS/bottle/bval_bottle.txt 

**Portal, Arizona**

S. K. Morgan Ernest, Glenda M. Yenni, Ginger Allington, Ellen K. Bledsoe, Erica M. Christensen, Renata M. Diaz, Keith Geluso, Jacob R. Goheen, Qinfeng Guo, Edward Heske, Douglas Kelt, Joan M. Meiners, Jim Munger, Carla Restrepo, Douglas A. Samson, Michele R. Schutzenhofer, Marian Skupski, Sarah R. Supp, Kate Thibault, Shawn Taylor, Ethan White, Hao Ye, Diane W. Davidson, James H. Brown, Thomas J. Valone. The Portal Project: a long-term study of a Chihuahuan desert ecosystem. bioRxiv 332783; doi: https://doi.org/10.1101/332783 
