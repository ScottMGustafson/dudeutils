from configparser import ConfigParser


def set_config_defaults(dct):
    if not "n" in list(dct['config'].keys()):
        dct['config']["n"] = 2


def tf(val):
    if val.lower() == 'true':
        return True
    elif val.lower() == 'false':
        return False
    elif val.isnumeric():
        return float(val)
    else:
        return val


def _set_values(dct):
    Config.ab_cfg = dct
    Config.cont_cfg = Config.ab_cfg.pop('continuum', {})
    Config.glob = Config.ab_cfg.pop('config', {})
    Config.vdisp = float(Config.glob.get('vdisp', 2.14))
    Config.vsig = float(Config.glob.get('vsig', 3.4))


class Config(object):
    ab_cfg = {}
    cont_cfg = {}
    glob = None
    vdisp = 2.14
    vsig = 3.4
    config_file = None

    @staticmethod
    def configure(config_file):
        """
        parse the config file to produce a dict of absorbers, parameters and 
        allowed ranges
    
        Input:
        ------
        config_file: config file to use.  data/random_sampling_config.cfg by default
    
    
        Output:
        -------
        dict of parameters and allowed ranges
    
        Raises:
        -------
        None
    
        """
        Config.config_file = config_file
        config = ConfigParser()
        config.optionxform = str
        config.read(config_file)
        dct = {}
        for item in list(config.sections()):
            dct[item] = dict(config.items(item))

        for item in dct.keys():
            if item == 'config':
                for key in dct[item].keys():
                    dct[item][key] = tf(dct[item][key])

            if item == 'continuum':
                for key, val in dct[item].items():
                    val = val.split(',')
                    if len(val) == 2:
                        dct[item][key] = {'xlim': float(val[0]), 'ylim': float(val[1])}
                    dct[item][key] = {'ylim': float(val[0])}
            else:
                for k in ["N", "b", "z"]:
                    try:
                        dct[item][k] = dct[item][k].replace(" ", "")
                        dct[item][k] = list(map(float, dct[item][k].strip().split(',')))
                    except KeyError:
                        pass

                cond = []
                if "N" in dct[item].keys():
                    cond += [dct[item]["N"][0] < 10., dct[item]["N"][-1] > 23.,
                             dct[item]["N"][-1] < dct[item]["N"][0]]

                if "b" in dct[item].keys():
                    cond += [dct[item]["b"][0] < 0., dct[item]["b"][-1] < dct[item]["b"][0]]

                if "z" in dct[item].keys():
                    cond += [dct[item]["z"][-1] < dct[item]["z"][0]]

                if any(cond):
                    raise Exception("check your random sampling inputs")

        set_config_defaults(dct)
        _set_values(dct)
