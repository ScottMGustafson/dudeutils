def get_keys(root, tag):
    children = root.find('CompositeSpectrum').findall(tag)
    if children==[]:
        raise Exception(str(children)+"does not contain "+str(tag))
    return list(dict(children[0].attrib).keys())

def prettify(elem):
    return et.tostring(elem, encoding="unicode")

def get_node(parent,tag,id):
    results=[item for item in parent.findall(tag) if item.get('id')==id]
    if len(results)>1:
        warnings.warn("more than one element contains the id \'%s\'\n returning the first one as default"%id)
        return results[0]
    elif len(results)==1:
        return results[0]
    else:
        raise Exception('id '+id+' not found')



