# Mitchell September 2017
# This script calculates the highest value from an arcgis tabulate area and creates a new table of proper format
# for input into Fleming et all script.
import arcpy, os
lyr =  'D:/ebird/Canada/BCR_Tabular.dbf'
fields = ['BCR_3', 'BCR_4', 'BCR_5', 'BCR_6', 'BCR_7', 'BCR_8', 'BCR_9', 'BCR_10', 'BCR_11', 'BCR_12', 'BCR_13', 'BCR_14', 'BCR_23', 'BCR_100']
newrow = {}
if not os.path.isfile('D:/ebird/Canada/CABCR_Table.dbf'):
  newtbl = arcpy.CreateTable_management('D:/ebird/Canada/' , 'CABCR_Table.dbf') 
  arcpy.AddField_management(newtbl, 'fips', "TEXT", field_length = 200)
  arcpy.AddField_management(newtbl, 'BCR', "TEXT", field_length = 200)

i = 0
for row in arcpy.SearchCursor(lyr):
  maxrow = 0
  maxbrc = ''
  for f in fields:
    if int(row.getValue(f)) > maxrow:
      maxrow = row.getValue(f)
      maxbcr = f[4:]
  newrow[i] = [row.getValue('DEGBLK'),maxbcr]
  i +=1

inscur = arcpy.InsertCursor("D:/ebird/Canada/CABCR_Table.dbf")
for k, v in newrow.iteritems():
  inrow = inscur.newRow()
  inrow.setValue('fips', v[0])
  inrow.setValue('BCR', v[1])
  inscur.insertRow(inrow)

del inscur
del inrow
