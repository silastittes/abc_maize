
for i in atypical_abc_out/*Na__1000_*; do grep -h -B2 "SFS:" $i | grep "SFS:" | cut -d ":" -f2- | sed 's/^ //g'; done > atypical_pop_sfs.txt &

for i in atypical_abc_out/*Na__1000_*; do grep -h -B2 "SFS:" $i | grep -h -A1 "#ngenes" | grep -v "\-\-" | grep -v "#Na"; done > atypical_pop.txt


#for i in `ls abc_out/*Na__200_*`; do grep -h -B2 "SFS:" $i | grep "SFS:" | cut -d ":" -f2- | sed 's/^ //g'; done > pop200_sfs.txt &

#for i in `ls abc_out/*Na__200_*`; do grep -h -B2 "SFS:" $i | grep -h -A1 "#Na" | grep -v "\-\-" | grep -v "#Na"; done > pop200.txt
