import dudeutils
import xml.etree.ElementTree as et


"""check that the format of an xml works.  if it doesn't correct and check for key collisions"""


    
pool = dict()


if __name__ == "__main__":
    tree = et.parse(fname)
    root = tree.getroot()
    for item in root:
        if "Lists" in item.tag:
            for it in item:
                try:
                    assert(it.get("id") not in pool.keys)
                    pool[it.get("id")]=it  #append node element to pool
                except:
                    #test for degeneracy in pool values.  if yes, simply delete node.
                    if not it in pool.values():
                        raise KeyError("critical error: key collision")
                    et.remove(it)
    tree.write(fname)
            

    
