# Scientific data provenance

`raw/Liefer2019_PANGAEA900248.zip` and its extracted Excel workbook are from
https://store.pangaea.de/Publications/Liefer-etal_2019/Final_MASS_Database.zip
Downloaded2026-10-16, CC-BY4.0. Cite:
Liefer, J.D.; Garg,A.; Fyfe,M.H.; Irwin,A.J.; Benner,I.; Brown,C.M.;
Follows,M.J.; Omta,A.W.; Finkel,Z.V. (2026). Elemental and macromolecular
content of marine diatoms and prasinophytes. PANGAEA.
https://doi.org/10.1594/PANGAEA.900248

Associated article: Liefer etal2019. The Macromolecular Basis of Phytoplankton
C:N:P Under Nitrogen Starvation. Front.Microbiol.10:763.
https://doi.org/10.3389/fmicb.2019.00763

Archive ZIP SHA256:
a91320274acff049b1acd71040ff6a0f2f9a6a8ef9af11e06397e23673b0782a

The landing page explicitly states that communication with the authors ceased
before usual final formal approval. Use the curation caveat. Six late-stage
replicates (O.tauri and T.pseudonana) are excluded for bacterial contamination
identified in the article. Blank=unmeasured; zero=true zero. The DIN unit label
in the workbook contradicts itself; article methods establish micromolar, equal
to mmolN/m3. O.tauri day2 C and N quotas duplicate day0 exactly. Dates are not
used to order experiments; days since N-free dilution are used instead.

`processed/culture_replicates.csv` contains81 rows,75 after the contamination
exclusion. It is extracted without altering the source workbook.
`processed/cultures.json` contains timepoint means/SD, source hash, protocol,
calibration/evaluation split, and a check of16 numerical Table1 values.
Maximum difference to rounded published values is0.396percent.

Initial post-dilution sampling day: Tp1, Ot2, Tw3, Micromonas6. The latter
starts at DIN0; it tests starvation persistence, not nutrient uptake kinetics.
Day0 replete C,N and measured cell volumes are used as reference traits even
for withheld species. The unknown earlier dilution is not guessed.
