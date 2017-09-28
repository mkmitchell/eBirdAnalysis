# Mitchell September 2017
# This script calculates the highest value from an arcgis tabulate area and creates a new table of proper format
# for input into Fleming et all script.
import arcpy, os
lyr =  'D:/ebird/Canada/JV_Tabular.dbf'
fields = ["EASTERN_HA", "WESTERN_BO", "PRAIRIE_PO", "UPPER_MISS", "PRAIRIE_HA", "CANADIAN_I", "ATLANTIC_C", "PACIFIC_BI", "INTERMOUNT"]
newrow = {}
if not os.path.isfile('D:/ebird/Canada/CAJV_Table.dbf'):
  newtbl = arcpy.CreateTable_management('D:/ebird/Canada/' , 'CAJV_Table.dbf') 
  arcpy.AddField_management(newtbl, 'fips', "TEXT", field_length = 200)
  arcpy.AddField_management(newtbl, 'JV', "TEXT", field_length = 200)
  arcpy.AddField_management(newtbl, 'AreasqKm', "DOUBLE")
  arcpy.AddField_management(newtbl, 'PartsqKm', "DOUBLE")
  arcpy.AddField_management(newtbl, 'propsqkm', "DOUBLE")

i = 0
for row in arcpy.SearchCursor(lyr):
  for f in fields:
    if int(row.getValue(f) > 1.0):
      print 'yes', f
      newrow[i] = [row.getValue('DEGBLK'),f, row.getValue('Total'), row.getValue(f), row.getValue(f)/row.getValue('Total')]
      i +=1

inscur = arcpy.InsertCursor("D:/ebird/Canada/CAJV_Table.dbf")
for k, v in newrow.iteritems():
  inrow = inscur.newRow()
  inrow.setValue('fips', v[0])
  inrow.setValue('JV', v[1])
  inrow.setValue('AreasqKm', v[2])
  inrow.setValue('PartsqKm', v[3])
  inrow.setValue('propsqkm', v[4])
  inscur.insertRow(inrow)

del inscur
del inrow
