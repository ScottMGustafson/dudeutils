import os
from shutil import copyfile

from dudeutils.data_types import Absorber, ContinuumPoint


class ModelIO(object):
    def __init__(self, path, output_name=None, overwrite=True):
        if not output_name:
            output_name = "output.csv"
        self.output = os.path.join(path, output_name)
        self.path = path
        if os.path.isfile(self.output):
            if overwrite:
                copyfile(self.output, self.output.replace('.csv', '_copy.csv'))
                os.remove(self.output)
            else:
                self.unique_name()

    def rename_output(self, n):
        path, fname = os.path.split(self.output)
        fname = "%d_output.csv" % n
        self.output = os.path.join(path, fname)

    def unique_name(self):
        i = 0
        self.rename_output(i)
        while os.path.isfile(self.output):
            self.rename_output(i)

    @staticmethod
    def ab_column_names(lst):
        cols = [(ab.id + "_ionName", ab.id + "_N", ab.id + "_b", ab.id + "_z") for ab in lst]
        return [x for sublst in cols for x in sublst]

    @staticmethod
    def ab_data(lst):
        return ["%s,%10.9lf,%10.9lf,%16.15lf" % (ab.ionName, ab.N, ab.b, ab.z) for ab in lst]

    @staticmethod
    def cnt_column_names(lst):
        cols = [(cnt.id + "_x", cnt.id + "_y") for cnt in lst if cnt.id != 'null']
        return [x for sublst in cols for x in sublst]

    @staticmethod
    def cnt_data(lst):
        return ["%15.11lf,%10E" % (cnt.x, cnt.y) for cnt in lst if cnt.id != 'null']

    @staticmethod
    def header_format(ab_lst, cnt_lst):
        abs = ModelIO.ab_column_names(ab_lst)
        conts = ModelIO.cnt_column_names(cnt_lst)
        return ','.join(['chi2'] + conts + abs) + '\n'

    def write(self, samples):
        self.header = ModelIO.header_format(samples[0][0], samples[0][1])
        inp = []
        for i in range(len(samples)):
            inp.append(','.join(
                ["%10.8lf" % (samples[i][2])] + ModelIO.cnt_data(samples[i][1]) + ModelIO.ab_data(
                    samples[i][0])) + '\n')

        if os.path.isfile(self.output):
            with open(self.output, 'a') as f:
                f.writelines(list(set(inp)))
        else:
            with open(self.output, 'w') as f:
                f.write(self.header)
                f.writelines(list(set(inp)))

    @staticmethod
    def _join_csv(path, prefix, output=None):
        if not output:
            output = os.path.join(path, prefix + '.csv')
        all_data = []
        header = ""
        for item in os.listdir(path):
            if os.path.isfile(os.path.join(path, item)) and item.startswith(prefix):
                with open(os.path.join(path, item), 'r') as f:
                    data = f.readlines()
                    if len(data) <= 1:
                        os.remove(os.path.join(path, item))
                        continue
                    header = data.pop(0)
                    all_data.extend(list(set(data)))
                os.remove(os.path.join(path, item))

        with open(output, 'w') as f:
            if len(header) > 0:
                f.write(header)
            else:
                raise Exception('too many headers.  check %s' % output)
            for line in sorted(list(set(all_data))):
                f.write(line)

    def join_csv(self, prefix='output_'):
        ModelIO._join_csv(self.path, prefix, output=self.output)

    @staticmethod
    def join_aborted_runs(path, prefix=None):
        def strip_prefix(fname):
            fname = fname.replace('.csv', '')
            parts = fname.split('_')
            if parts[-1].isnumeric() and parts[-2].isnumeric():
                return '_'.join(parts[:-2])
            else:
                return fname

        def get_prefix_set():
            return list(set([strip_prefix(it) for it in list(os.listdir(path)) if it.endswith('.csv')]))

        if prefix is None:
            lst = get_prefix_set()
        else:
            lst = [prefix]
        for pfix in lst:
            try:
                ModelIO._join_csv(path, pfix)
            except:
                print(pfix, path)
                raise

    @staticmethod
    def find_best(lines: list) -> list:
        """
        
        Parameters
        ----------
        lines : list
            list of strings, unparsed csv data

        Returns
        -------
        list
            best model in format of [chi2, cont_pts , .. , absorbers, .. ] 
        """
        min_chi2 = None
        curr = None
        for line in lines[1:]:
            lst = line.strip().split(',')
            if not min_chi2:
                try:
                    min_chi2 = float(lst[0].strip())
                except:
                    print(lst)
                    raise
                curr = lst
            elif float(lst[0]) < min_chi2:
                min_chi2 = float(lst[0].strip())
                curr = lst
        return curr

    @staticmethod
    def lst_to_model(model, header, lst):
        """
        
        Parameters
        ----------
        model : model.Model
            model should have the same skeleton as the desired output, i.e. same number of absorbers and continuum 
            points as lst
        header : list
            a list of header keys from the input data
        lst : list
            model to parse from csv data

        Returns
        -------
        model.Model 
            Model instance parsed from lst
        """
        model = model.copy()
        dct = dict(zip(header, lst))

        cont_lst, ab_lst = [], []
        model.chi2 = float(lst[0])
        ids = list(set([k.split('_')[0] for k in header[1:]]))
        for _id in ids:
            subdct = {k.split('_')[-1]: val.strip() for k, val in dct.items() if k.startswith(_id)}
            if 'ionName' in subdct.keys():
                ab_lst.append(Absorber(**subdct))
            elif 'x' in subdct.keys():
                cont_lst.append(ContinuumPoint(**subdct))

        model.set_absorbers(ab_lst)
        model.set_cont_points(cont_lst)
        return model

    @staticmethod
    def read_best(model, fname):
        with open(fname, 'r') as f:
            assert os.stat(fname).st_size > 0, 'Empty File'
            lines = f.readlines()

        header = lines.pop(0).strip().split(',')
        best = ModelIO.find_best(lines)
        return ModelIO.lst_to_model(model, header, best)

# if __name__ == '__main__':
#     path='/home/scott/J0744/'
#     ModelIO(path=path).join_csv(prefix)
