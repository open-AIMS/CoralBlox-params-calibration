# Historical disturbance (storms/COTS) mortality rates

- **Dataset custodian** : Pedro Ribeiro de Almeida
- **Contact** : p.ribeirodealmeida@aims.gov.
- **Purpose**: This dataset is for academic research purposes only. This dataset was generated as part of the CoralBlox-ADRIA calibration process.
- **Date published** :  2025 11 05

---


# Data structure

This dataset contains a single NetCDF file. The variable `disturbance_mortality_scens`
has dimensions:

- timesteps : Vector{Int64} 2008:2022,
- locs : Vector{String} ['10-330', …, '23-049'],
- species : Vector{String} ['tabular_Acropora', 'corymbose_Acropora', 'corymbose_non_Acropora',
'small_massives', 'large_massives'],
- scenarios : Vector{Int64} [1,2] (both scenarios are identical),

Values represent proportion of coral cover lost (refered to here as "mortality rate") (0-1)
per year, reef (loc), and functional group (species) after a storm or COTS event. Both
scenarios are identical. The mortality rates for each functional group are relative to each
species cover thus they they don't sum up to one for a particular year and reef.


# Description

This dataset represents mortality rates after a storm or COTS event derived from historical
coral cover data observed by LTMP per functional group, year and each reef of the GBR. Years
range from 2008 to 2022; functional groups represented are tabular Acropora, corymbose
Acropora, corymbose non Acropora, small massives, large massives. Scenarios 1 and 2 are
identical (this is a limitation of the tool used to create the NetCDF files).

Disturbance and coral cover data used to generate this dataset were extracted from the Reef
Monitoring Dashboard. The steps followed to create this dataset were:

- Generate a list of storm and/or COTS disturbances within the considered timeframe for
each reef with available data. The process of determining the years before and after was
done manually.
- For each disturbance in each reef, take the relative coral cover before and after that,
and use that to calculate an overall mortality rate (accounting for all functional groups
together)
- For each disturbance in each reef, take the relative coral composition before and after
that and use that, together with the overall mortality rate, to calculate the mortality rate
for each functional group
- For reefs with no coral composition data, compare the average change in coral composition
for across the whole GBR with the average change in coral composition across the reef on
the same bioregion and use the one with the lower coefficient of variation (this happened
for 24 of the 154 pairs disturbance/reef)

## Caveats

The dimension `timesteps` represents years when a mortality rate should be applied. However,
in some cases there is a gap of more than one year between the years with observed cover
before and after a disturbance. In those cases, we gave preference to aligning the data on
the year after the disturbance.

For two reefs ("Rebe Reef" and "Martin Reef") there are disturbances reported in the
transect data (coral composition) that are not reported in the manta tow data, causing a big
time window in the manta tow without data reported. In Rebe Reef there are two manta tow
disturbances, one in 2005 and another in 2011, but in the transect data for the same reef
there is only a disturbance in 2009. The same happens in "Martin Reef" between 2013 and 2017,
with an extra disturbance reported in the transect data in 2015. For these, a new datapoint
was manually added to have a "hook" where the data for this extra disturbance could be added.
The same happened on the other way for "Taylor Reef" between 2016 and 2018 and for
"Hoskyn Island" between 2008 and 2010, where there is a disturbance reported in the manta
tow data that was not reported in the transect. In those cases an extra year was manually
added to the transect dataset and the mortality was "split" in two years.

We assume that the maximum mortality rate is 0.95. That is to handle cases where a
functional group cover recorded data goes from a certain positive value to 0 - which can be
misleading, since a very low density of some functional group could lead to no cover being
recorded for that functional group just by chance.

As mentioned, missing benthic disturbance cover data for some specific combinations of
reef/disturbance is filled with either average across all reefs or within bioregion. That
(missing benthic disturbance cover data) happened on 24 of the 154 pairs reef/disturbance.

# License
[Copyright - All rights reserved](https://docs.provena.io/licenses.html#copyright-all-rights-reserved)
