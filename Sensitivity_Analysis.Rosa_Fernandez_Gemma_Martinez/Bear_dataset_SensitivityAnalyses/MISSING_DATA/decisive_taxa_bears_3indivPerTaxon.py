import fileinput, glob, string, sys, os, csv
import shutil

#it counts the occurrences of species names, sums these for taxa and then prints to the screen those orthogroups that are present only in certain predefined taxa or groups of taxa; the respective decisive OGs need then to be copied to a new folder

orthogroup = glob.glob("*.fa")   #change extension accordingly
dirname_Decisive = 'Decisive_genes3'
dirname_NonDecisive = 'NonDecisive_genes3'

os.mkdir(dirname_Decisive)
os.mkdir(dirname_NonDecisive)

for filename in orthogroup:
    fh = open(filename)
    content = fh.read()
    fh.close()
    
    Maria_count = content.count("Maria")
    Noah_count = content.count("Noah")
    Margo_count = content.count("Margo")
    Summer_count = content.count("Summer")
    Siro_count = content.count("Siro")
    Montana_count = content.count("Montana")
    Juan_count = content.count("Juan")
    Adelaide_count = content.count("Adelaide")
    Pepe_count = content.count("Pepe")
    Amparo_count = content.count("Amparo")
    Joseph_count = content.count("Joseph") 
    Luisa_count = content.count("Luisa") 
    Margaret_count = content.count("Margaret")
    Maripepa_count = content.count("Maripepa")
    Paco_count = content.count("Paco")
    Oskar_count = content.count("Oskar")
    
 
    Ailuropoda_sum = Luisa_count + Pepe_count + Juan_count + Siro_count
    UrsusMaritimus_sum = Maria_count + Maripepa_count + Margaret_count + Joseph_count
    UrsusArctos_sum = Margo_count + Paco_count + Adelaide_count + Amparo_count
    UrsusAmericanus_sum = Noah_count + Montana_count + Summer_count + Oskar_count
        
# in the following groups of taxa are created that contain each gene at least once each, and the gene should be misisng in all other groups; results are to be printed to screen
    if Ailuropoda_sum >= 3 and UrsusMaritimus_sum >= 3 and UrsusArctos_sum >= 3 and UrsusAmericanus_sum >= 3:
        print("Decisive", filename)
        shutil.copy(filename, dirname_Decisive)
	
    if Ailuropoda_sum < 3 or UrsusMaritimus_sum < 3 or UrsusArctos_sum < 3 or UrsusAmericanus_sum < 3:
        print("Not_Decisive", filename)
        shutil.copy(filename, dirname_NonDecisive)
	
	  

 
