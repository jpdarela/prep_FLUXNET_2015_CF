from nc_write import write_ref_data, write_site_nc, create_gridlist
from conversions import mod_pr, mod_tas

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
    
    # create modified PREC
    write_site_nc('pr', mod=mod_pr)
    write_site_nc('tas', mod=mod_tas)
    
if __name__ == "__main__":
    main()
