import sys
import xml.etree.ElementTree as et
import dude_xmlutils
import astronomy_utils as astro

""" 
    command line script to write new values to xml or to dump current vals into 
    text file.
    a typical command will be like:
      `> python cmd_line_io.py get file.xml id=foobar id2=foo id3=bar`
    this will dump all the data from the listed absorbers to a file 
    `absorberData.txt`.

      `> python cmd_line_io.py write file.xml id=foobar N=12.5 b=8.3 z=4.0`
    will right the fitting parameters back into the xml
    

"""

def parser(arg_list):
  datalst = []
  i=0
  action = 'write' if 'write' in arg_list else 'get'
  while i < len(arg_list):
    if '.xml' in arg_list[i]:
      xml_file = arg_list[i]
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
        datalst.append(dude_xmlutils.Data(**kwargs))
      elif key=='chisq':
        chisq=int(arg_list[i].split('=')[1])
    i+=1

  return action, xml_file, datalst, chisq

action, xml_file, datalst, chisq = parser(sys.argv[1:])

for i in range(1,len(datalst)):
  setattr(datalst[i],'vel',astro.get_vel_shift(datalst[0].z,datalst[i].z)[0])
  
if action=='get': 
  for item in datalst:
    item.getData()
  output_string = ''
  for item in datalst:
    output_string+=str(item)
  output_string+=' chisq='+str(chisq)
  f = open('absorberData.txt','a') # append read data to file
  f.write(output_string+'\n')
  f.close()

elif action=='write':
  for item in datalst:
    item.writeData()
  tree.write(xml_file)

else:
  raise Exception('unrecognized command: '+str(action))

