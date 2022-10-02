from nc_write import write_ref_data, write_site_nc, create_gridlist

# Forcing data
VAR = 'tas', 'vpd', 'hurs', 'ps', 'pr', 'wind', 'rsds'

# Reference variables
REF = 'gpp', 'reco', 'nee', 'aet'

# FLUXNET 2015 Sites
SITES =  ['Dav', 'Tha', 'Hai', 'Lnf', 
          'Obe', 'Lae', 'BK1', 'Lkb', 
          'Sor', 'Col', 'Ren', 'Fyo']

def main():
    for variable in VAR:
        write_site_nc(variable)
        
    create_gridlist("FLUXNET2015_gridlist")
    
    for variable in REF:
        for site in SITES:
            write_ref_data(variable, site)

if __name__ == "__main__":
    main()
