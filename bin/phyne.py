
import os
import json

class PHYNE:
    def __init__(self, args):
        self.mode = args.mode
        self.config = args.config
        self.outdir = args.outdir
        self.prefix = args.prefix

    def load_conf_json(self):
        self.conf_dic = json.load(open(self.config))

    def make_cmd_mlst_with_scheme(self):
        cmd = ['python']
        cmd.append(self.conf_dic["phyne_mlst_exe"])
        cmd.append('--target-scheme')
        cmd.append(self.conf_dic["target_scheme"])
        label_s = list()
        infn_s = list()
        for label, infn in self.conf_dic["input_genome"].items():
            label_s.append(label)
            infn_s.append(infn)
        cmd.append('--infn-s')
        cmd.extend(infn_s)
        cmd.append('--label-s')
        cmd.extend(label_s)
        cmd.append('--outdir')
        cmd.append(self.outdir)
        cmd.append('--prefix')
        cmd.append(self.prefix)
        return cmd

    def make_cmd_mlst_without_scheme(self):
        cmd = ['python']
        cmd.append(self.conf_dic["phyne_mlst_exe"])
        label_s = list()
        infn_s = list()
        for label, infn in self.conf_dic["input_genome"].items():
            label_s.append(label)
            infn_s.append(infn)
        cmd.append('--infn-s')
        cmd.extend(infn_s)
        cmd.append('--label-s')
        cmd.extend(label_s)
        cmd.append('--outdir')
        cmd.append(args.outdir)
        cmd.append('--prefix')
        cmd.append(args.prefix)
        return cmd

def main(args):
    phyne = PHYNE(args)
    if args.mode in ['mlst']:
        phyne.load_conf_json()
        if phyne.conf_dic["target_scheme"]:
            cmd = phyne.make_cmd_mlst_with_scheme()
        else:
            cmd = phyne.make_cmd_mlst_without_scheme()
        print (' '.join(cmd))
        os.system(' '.join(cmd))

    elif args.mode in ['nssnp']:
        pass

    elif args.mode in ['ortholog']:
        pass



if __name__=='__main__':
    import argparse
    parser = argparse.ArgumentParser()
    parser.add_argument('--mode', choices=('mlst', 'nssnp', 'ortholog'))
    parser.add_argument('--config')
    parser.add_argument('--outdir')
    parser.add_argument('--prefix')
    args = parser.parse_args()
    main(args)
