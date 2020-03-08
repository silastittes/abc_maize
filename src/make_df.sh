grep -h -B2 "SFS:" abc_out/* | grep "SFS:" | cut -d ":" -f2- | sed 's/^ //g' > pop_sfs.txt

grep -h -B2 "SFS:" abc_out/* | grep -h -A1 "#Na" | grep -v "\-\-" | grep -v "#Na" > pop.txt
