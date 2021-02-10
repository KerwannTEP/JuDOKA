using HDF5 # To have access to .hf5 files
using Plots # To be able to plot data

########################################
# Function that read the .hf5 file
########################################

namefile = "../../data/Dump_Diffusion_Coefficients.hf5"

"""
    openData(namefile)

Recovers the data from the .hf5 file `namefile` and returns it as a 4-uple.
In the order nbj, aMeasure, tabj, tabDRRjj.
"""
function openData(namefile)
    file = h5open(namefile, "r")
    nbj = read(file,"nbj")
    nba = read(file,"nba")
    tabaj = read(file,"tabaj")
    tabj = read(file,"tabj")
    taba = read(file,"taba")
    tabDRRjj = read(file,"tabDRRjj")
    close(file)
    return nbj, nba, tabaj, tabj, taba, tabDRRjj
end

########################################
# Getting the data
########################################
println("Recovering plot data...")
nbj, nba, tabaj, tabj, taba, tabDRRjj = openData(namefile)

########################################
# Converting the data into log-scaling
########################################
println("Converting data into log-log scaling...")
for ij=1:nbj
    tabj[ij] = log10(tabj[ij])
end

for ia=1:nba
    taba[ia] = log10(taba[ia])
end


for iGrid=1:nba*nbj
    if (tabDRRjj[iGrid] == 0) # Within loss-cone, set an arbitrary low value for null
        tabDRRjj[iGrid] = -10
    else
        tabDRRjj[iGrid] = log10(tabDRRjj[iGrid])
    end
end

########################################
# Plotting the data
########################################
println("Plotting the data...")
p = Plots.contourf(tabj,taba,tabDRRjj)
Plots.savefig(p,"../../graphs/Julia/Djj.png") # Saves the figure
Plots.display(p)
readline()
