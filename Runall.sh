python compare_cupSODA_scipy_lsoda.py tyson new_tyson_results.log > new_tyson_log.txt
python parse_cupsoda_log.py new_tyson_log.txt new_tyson_results.log TYSON.txt
python compare_cupSODA_scipy_lsoda.py ras new_ras_results.log > new_ras_log.txt
python parse_cupsoda_log.py new_ras_log.txt new_tyson_results.log RAS.txt
python compare_cupSODA_scipy_lsoda.py earm new_earm_results.log> new_earm_log.txt
python parse_cupsoda_log.py new_earm_log.txt new_earm_results.log EARM.txt
cat TYSON.txt > GPU_timing.csv
cat RAS.txt >> GPU_timing.csv
cat EARM.txt >> GPU_timing.csv

