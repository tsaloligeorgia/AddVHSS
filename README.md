## VHSS  

The instruction is applied for Arch Linux. 

### REQUIRED LIBS

gmp gmpxx 

### INSTALLATION 

pacman -S gmp 

### BUILD

Select the desired version of DAGS you would like to build by editing
the makefile 
To build the binary, just run:

```bash
cd add_vhss_hss
make clean 
make 
cd ..

cd add_vhss_hss
make clean 
make 
cd ..

cd add_vhss_hss
make clean 
make 
cd .. 

```


### RUN 
Select which one you desire to run. 
Go to add_vhss_hss to VAHSS-HSS
```bash
cd add_vhss_hss
./add_vhss_hss
```
Go to add_vhss_lhs to VAHSS-LHS
```bash
cd add_vhss_lhs
./add_vhss_lhs
```
Go to add_vhss_tss to VAHSS-TSS
```bash
cd add_vhss_tss
./add_vhss_tss
```
