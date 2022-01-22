ptmin=100
ptmin=80

./leadsj_vs_zloss.py --jetptmin ${ptmin} --jetptmax 120 --nev 2000 --py-hardQCD --out lsj_zloss_any.root
./leadsj_vs_zloss.py --jetptmin ${ptmin} --jetptmax 120 --nev 2000 --py-hardQCD --out lsj_zloss_any_kt.root --kt

./leadsj_vs_zloss.py --jetptmin ${ptmin} --jetptmax 120 --nev 2000 --py-hardQCDgluons --out lsj_zloss_glue.root
./leadsj_vs_zloss.py --jetptmin ${ptmin} --jetptmax 120 --nev 2000 --py-hardQCDgluons --out lsj_zloss_glue_kt.root --kt

./leadsj_vs_zloss.py --jetptmin ${ptmin} --jetptmax 120 --nev 2000 --py-hardQCDquarks --out lsj_zloss_quark.root
./leadsj_vs_zloss.py --jetptmin ${ptmin} --jetptmax 120 --nev 2000 --py-hardQCDquarks --out lsj_zloss_quark_kt.root --kt
