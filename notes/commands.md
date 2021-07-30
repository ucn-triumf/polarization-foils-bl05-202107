# compress the data files
tar -czvf data_arch/data_210721.tar.gz data/

# extract the data 
# place data_210721.tar.gz to your repository 
tar -xvzf data_210721.tar.gz 


# ROOT commands 
- open root files 
  >  root -l data/210713_SiFe/20210717055945_list.root
- List the elements included in the .root file
  > .ls
- Print information of TTree named "T"
  > T->Print()
- Print
  > T->GetEntries()
