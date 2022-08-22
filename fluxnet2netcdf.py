from nc_write import write_site_nc, save_gridlist


VAR = 'tas', 'vpd', 'hurs', 'ps', 'pr', 'co2', 'wind', 'rsds'
SITE = 'Dav', 'Tha'




def main():
    for variable in VAR:
        for site in SITE:
            write_site_nc(variable, site)
    
    for site in SITE:
        save_gridlist(site)

if __name__ == "__main__":
    main()
    