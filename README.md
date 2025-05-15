# ACE 592 Final Project - Neighborhood Hotspot Analysis

This repository contains the code and data for the ACE 592 Final Project, which focuses on identifying and analyzing crime hotspots in Chicago neighborhoods using crime and ShotSpotter data. The project aims to map and analyze spatial patterns of violent crime, particularly homicides, in relation to demographic and socioeconomic data.

## Project Objective
The objective of this project is to investigate the spatial distribution of violent crime across Chicago neighborhoods, focusing on identifying high-risk areas ("hot spots") using spatial analysis techniques. The analysis explores the correlation between crime data and ShotSpotter alerts within a defined proximity range.

## Data Sources
- [Violence Reduction - Victims of Homicides and Non-Fatal Shootings](https://data.cityofchicago.org/Public-Safety/Violence-Reduction-Victims-of-Homicides-and-Non-Fa/gumc-mgzr/about_data)
- [Census Data - Selected Socioeconomic Indicators in Chicago](https://data.cityofchicago.org/Health-Human-Services/Census-Data-Selected-socioeconomic-indicators-in-C/kn9c-c2s2/about_data)
- [Violence Reduction - ShotSpotter Alerts Historical](https://data.cityofchicago.org/Public-Safety/Violence-Reduction-Shotspotter-Alerts-Historical/3h7q-7mdb/about_data)
- [Transportation Data](https://data.cityofchicago.org/Transportation/transportation/7ez8-272k/about_data)
- [Chicago Maps and Spatial Data](https://data.cityofchicago.org/browse?q=map&sortBy=relevance&pageSize=20&page=1)

## Key Analysis Steps
- Data Preprocessing: Cleaning and merging datasets, converting to GeoDataFrames.
- Spatial Analysis: Identifying hotspots using proximity analysis.
- Visualization: Mapping hotspots and creating heatmaps for visual representation.

## Outputs
- **Spatial Analysis Map:** Hotspots visualized on the Chicago map with neighborhood boundaries.
- **Correlation Analysis:** Statistical analysis of ShotSpotter alerts and homicide data.


## Contact
For any questions or further information, please contact Austin Belmonte at abelma2@illinois.edu.

## License
This project is licensed under the MIT License. See the LICENSE file for details.
