using HDF5 # To have access to .hf5 files
using Plots # To be able to plot data

########################################
# Function that read the .hf5 file
########################################

namefile = "../../data/Dump_Diffusion_Coefficients_Cut.hf5"

"""
    openCutData(namefile)

Recovers the data from the .hf5 file `namefile` and returns it as a 4-uple.
In the order nbj, aMeasure, tabj, tabDRRjj.
"""
function openCutData(namefile)
    file = h5open(namefile, "r")
    nbj = read(file,"nbj")
    aMeasure = read(file,"aMeasure")
    tabj = read(file,"tabj")
    tabDRRjj = read(file,"tabDRRjj")
    jlc = read(file,"jlc")
    close(file)
    return jlc, nbj, aMeasure, tabj, tabDRRjj
end

########################################
# Getting the data
########################################

println("Recovering plot data...")
jlc, nbj, aMeasure, tabj, tabDRRjj = openCutData(namefile)

########################################
# Removing diffusion coefficients within the loss-cone
########################################

ij = 1

while (tabj[ij] < jlc)
    global ij += 1
end

tabj = tabj[ij:nbj]
tabDRRjj = tabDRRjj[ij:nbj]

########################################
# Converting the data into log-scaling
########################################

println("Converting data into log-log scaling...")

for ij=1:length(tabj)
    tabj[ij] = log10(tabj[ij])
end

for ij=1:length(tabDRRjj)
    tabDRRjj[ij] = log10(tabDRRjj[ij])
end

########################################
# Plotting the data
########################################

println("Plotting the data...")

p = plot(tabj, tabDRRjj, legend=false)#, scale=:log10)
savefig(p,"../../graphs/Julia/DjjCut.png") # Saves the figure
display(p) # Display plot

readline() # Plot window stays open until we press "Enter"
