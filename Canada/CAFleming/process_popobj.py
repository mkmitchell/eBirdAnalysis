# Mitchell September 2017
# The purpose of this script is to convert Fleming et al. R script output (Method 4b abd 4d county pop80 and LTA totals to the orignal format.
# This output will be slightly modified to feed Mitchell R scripts for calculating population area under the curve by county and degree block.
####
import csv, os

workdir =  'D:/ebird/Canada/CAFleming'
for method in ['4B', '4D']:
  for poptype in ['LTA', 'pop80']:
    lyr = ''
    if poptype == 'pop80':
      lyr = os.path.join(workdir, 'Meth' + method + 'countypop80.csv')
    else:
      lyr = os.path.join(workdir, 'Meth' + method + 'countypopLTA.csv')
    species = ['ABDU', 'AGWT', 'AMWI', 'BUFF', 'BWTE', 'CANV', 'COEI', 'GADW', 'KIEI', 'LTDU', 'MALL', 'NOPI', 'NSHO', 'REDH', 'RNDU', 'RUDU', 'SCAU', 'WODU', 'CITE']
    popobj = {}
    i=0
    with open(lyr) as csvfile:
      mergelist = ['fips'] + species
      reader = csv.DictReader(csvfile, fieldnames=mergelist)
      for row in reader:
        if i == 0:
          i += 1
          continue
        for sp in species:
          popobj[i] = [row['fips'], sp, row[sp]]
          i += 1

    #write out
    with open(os.path.join(workdir, method + '_' + poptype + '.csv'), 'w') as outcsvfile:
      if poptype == 'pop80':
        popcol = '80percPopObj'
      else:
        popcol = 'LTAPopObj'
      fieldnames = ['fips', 'species', popcol, 'CODE']
      writer = csv.DictWriter(outcsvfile, fieldnames=fieldnames, delimiter=',', lineterminator='\n')
      writer.writeheader()
      for k, v in popobj.iteritems():
        writer.writerow({'fips': v[0], 'species': v[1],popcol: v[2], 'CODE': method})
