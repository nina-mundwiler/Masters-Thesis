## 🔧 Dependencies

The script requires the following R packages:

```r
install.packages(c("lubridate", "dplyr", "ggplot2", "zoo"))
```

## 📊 Data Access & Licensing

The data used in this project was obtained from the official portal of the Swiss Federal Office of Meteorology and Climatology (MétéoSwiss), specifically intended for **educational and research** purposes:

➡️ [MétéoSwiss – Data for Education and Research](https://www.meteosuisse.admin.ch/services-et-publications/service/produits-meteorologiques-et-climatiques/portail-de-donnees-pour-l-enseignement-et-la-recherche.html)

To access the data:

1. Visit the [IDAweb portal](https://gate.meteoswiss.ch/idaweb) (registration required)
2. Request data for NEU station (`NEU_rre150d0.txt`)
3. Comply with MétéoSwiss's usage policy

> 🔐 **Important:** The raw data file may not be redistributed or published in public repositories (e.g., GitHub) without explicit permission from MétéoSwiss. Only derived analyses and visualizations may be shared publicly.

