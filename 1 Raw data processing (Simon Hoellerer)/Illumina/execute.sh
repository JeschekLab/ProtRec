cd "/links/groups/panke/01_Exchange_Folder/02_PhD_students_and_seniors/40_hoellerer/Data_Management/Chapter_06_Lukas_TEV_NGS/uASPIre_screen4_5/"

./00_execute.sh 00_config_S4_1.cfg

./00_execute.sh 00_config_S4_2.cfg

./00_execute.sh 00_config_S4_3.cfg

Rscript ./combined/01_LH.R
Rscript ./combined/02_LH.R
