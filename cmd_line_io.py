import sys
import xml.etree.ElementTree as et
import xmlutils

def parser(arg_list):
  datalst = []
  i=0
  action = 'write' if 'write' in arg_list else 'get'
  while i < len(arg_list):
    if "=" in arg_list[i]:
      key, val = arg_list[i].split('=')
      if key=='id':
        kwargs = {}
        kwargs['id']=val
        i+=1
        while arg_list[i][:2]!='id' and arg_list[i][:5]!='chisq':
          key, val = arg_list[i].split('=')
          kwargs[key] = val
          i+=1
        datalst.append(Data(**kwargs))
      elif key=='chisq':
        chisq=int(arg_list[i].split('=')[1])
    if '.xml' in arg_list[i]:
      xml_file = arg_list[i]
    i+=1

  return action, xml_file, datalst, chisq

action, xml_file, datalst, chisq = parser(sys.argv[1:])
tree = et.parse(xml_file) #parse xml file
root = tree.getroot()

for i in range(1,len(datalst)):
  setattr(datalst[i],'vel',getVelocityShift(datalst[0].z,datalst[i].z)[0])
  
if action=='get': 
  for item in datalst:
    item.getData(root)
  output_string = ''
  for item in datalst:
    output_string+=str(item)
  output_string+=' chisq='+str(chisq)
  f = open('absorberData.txt','a') # append read data to file
  f.write(output_string+'\n')
  f.close()

elif action=='write':
  for item in datalst:
    item.writeData(root)
  tree.write(xml_file)

else:
  raise Exception('unrecognized command: '+str(action))

