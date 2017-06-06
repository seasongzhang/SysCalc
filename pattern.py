# -*- coding:utf-8 -*-

from __future__ import division
import numpy as np
import pandas as pd


class Pattern:
    def __init__(self, mtype, spd):
        """

        :param mtype: 曳引机型号，如"ZPML-J117","ZPML-GH70"
        :param spd: 额定速度，单位为m/min
        """

        self.spd = float(spd)

        # 根据曳引机的型号来判断是否是高速梯，并结合速度来选取圆角表格

        power = mtype.split('-')[1][1:4]

        try:
            is_low_spd = int(power) < 400
        except ValueError:
            is_low_spd = False
        if not is_low_spd:
            self.acc, self.dec = [0.8, 0.8]
            self.t = {'t1': 1.1, 't3': 1.0, 't5': 1.0, 't7': 1.2}
        elif spd >= 150:
            self.acc, self.dec = [0.8, 0.8]
            self.t = {'t1': 0.8, 't3': 0.8, 't5': 0.8, 't7': 0.8}
        else:
            self.acc, self.dec = [0.6, 0.6]
            self.t = {'t1': 0.6, 't3': 0.6, 't5': 0.6, 't7': 0.6}

    def calc_by_distance(self, distance):
        """根据在Mathmatica中的计算，针对唐启峰提供的3套加速度/减速度/圆角时间的参数，若仅保证4段圆角时间，3套参数对应的
        行程分别为0.432m/1.024m/1.736m。因此，只要保证层间距大于1.736m，单次行程就必然包含完整的4个圆角，这在一定程
        度上简化了模型（计算）。基于行程的计算流程为：
        a. 判断行程是否能保证4段圆角时间
        b. 判断计算所得的额定速度运行距离s_uniform是否大于0？
            b1. 若>0，直接计算匀速运行时间；
            b2. 若<0，说明在该行程下系统无法达到额定速度，得根据3phase的计算方式得到分段时间。
        
        :param distance: 
        :return: 
        """
        # 为计算公式方便，创建本地变量
        d = distance
        v = self.spd / 60
        acc = self.acc
        dec = self.dec
        t1 = self.t['t1']
        t3 = self.t['t3']
        t5 = self.t['t5']
        t7 = self.t['t7']
        t2 = t4 = t6 = 0

        if d < 1.736:  # todo: 针对不同的梯种设置不同的限制
            print(u'不允许单次行程小于1.736m！')
        else:
            t2 = v / acc - (t1 + t3) / 2
            t6 = v / dec - (t5 + t7) / 2
            s = [acc / 6 * t1 ** 2,  # t1内的行程
                 acc / 2 * (t1 + t2) * t2,  # t2内的行程
                 acc * (t1 / 2 + t2 + t3 / 3) * t3,  # t3内的行程
                 v * t5 - dec / 6 * t5 ** 2,  # t5内的行程
                 t6 * (v - dec / 2 * t5 - dec / 2 * t6),  # t6内的行程
                 dec / 6 * t7 ** 2  # t7内的行程
                 ]
            s_uniform = d - sum(s)
            if s_uniform >= 0:
                # s.insert(index=3, object=s_uniform)
                t4 = s_uniform / v
            else:
                t4 = 0
                delta = np.sqrt(72 * acc * d - 24 * dec * d - 3 * acc ** 2 * t1 ** 2
                                + acc * dec * t1 ** 2 + 6 * acc ** 2 * t3 ** 2
                                - acc * dec * t3 ** 2 + 6 * acc ** 2 * t3 * t5
                                + 3 * acc ** 2 * t5 ** 2 + 3 * acc * dec * t5 ** 2
                                - dec ** 2 * t5 ** 2 - 6 * acc ** 2 * t3 * t7
                                + 6 * acc * dec * t3 * t7 - 6 * acc ** 2 * t5 * t7
                                + 6 * acc * dec * t5 * t7 + 3 * acc ** 2 * t7 ** 2
                                - 9 * acc * dec * t7 ** 2 + 4 * dec ** 2 * t7 ** 2
                                )
                t2 = 0.166667 / (3 * acc - dec) * (-9 * acc * t1 + 3 * dec * t1
                                                   - 12 * acc * t3 + 3 * dec * t3
                                                   - 3 * acc * t5 + 3 * acc * t7
                                                   - 3 * dec * t7 + 1.73205 * delta)
                t6 = 0.5 / (3 * acc - dec) * (-acc * t3 - 3 * acc * t5
                                              - 2 * acc * t7 + 0.57735 * delta)

        # 计算后输出结果
        self.t['t4'] = t4
        self.t['t2'] = t2
        self.t['t6'] = t6
        self.t['t8'] = 8
        # print(delta, t2, t4, t6)
        periods = [self.t[k] for k in sorted(self.t.keys())]

        return [self.acc, self.dec] + periods + [self.calc_whole_distance()]

    def calc_by_tr(self, tr, load_ratio=0.5):
        """

        :param tr: 提升高度
        :param load_ratio: 负载持续率，根据提升高度计算得到电梯单次运行时间后根据负载持续率计算休息时间
        :return:
        """
        tr = tr
        self.calc_by_distance(distance=tr)
        t_total = sum(self.t.values()) / load_ratio
        self.t['t8'] = t_total * (1 - load_ratio)
        periods = [self.t[k] for k in sorted(self.t.keys())]
        return [self.acc, self.dec] + periods + [self.calc_whole_distance()]

    def calc_by_sph(self, sph=180, load_ratio=0.5):
        """
        由sph和load_ratio可以决定工作制
        :param sph: start per hour
        :param load_ratio: 负载持续率，无论sph值如何，保证每小时内有load_ratio%的时间电梯处于非停止状态
        :return:
        """
        v = self.spd / 60
        spd = self.spd
        acc = self.acc
        dec = self.dec
        t1 = self.t['t1']
        t3 = self.t['t3']
        t5 = self.t['t5']
        t7 = self.t['t7']

        t2 = spd / 60 / acc - (t1 + t3) / 2
        t6 = spd / 60 / dec - (t5 + t7) / 2

        t4 = load_ratio * 3600 / sph - sum(self.t.values())

        # 若计算得到的匀速时间>0，说明按照sph和负载持续率计算确实存在匀速运行阶段
        # 若计算得到的匀速时间<0，说明按照sph和负载持续率计算不存在匀速运行阶段，需要考虑加速过程中直接减速的情况
        if t4 >= 0:
            t4 = t4
        else:
            t4 = 0
            t2_plus_t6 = 3600 / sph - (t1 + t3 + t5 + t7)
            t2_minus_t6 = ((t5 + t7) - (t1 + t3)) / 2
            t2 = (t2_plus_t6 + t2_minus_t6) / 2
            t6 = (t2_plus_t6 - t2_minus_t6) / 2
        t8 = 3600 / sph * (1 - load_ratio)

        self.t['t4'] = t4
        self.t['t2'] = t2
        self.t['t6'] = t6
        self.t['t8'] = t8
        periods = [self.t[k] for k in sorted(self.t.keys())]
        return [self.acc, self.dec] + periods + [self.calc_whole_distance()]

    def calc_whole_distance(self):
        # 基于加减速度和各阶段时长，计算总行程
        acc = self.acc
        dec = self.dec

        wd = 0
        wd += acc / 6 * (self.t['t1'] ** 2)
        wd += acc / 2 * (self.t['t1'] + self.t['t2']) * self.t['t2']
        wd += acc * (self.t['t1'] / 2 + self.t['t2'] + self.t['t3'] / 3) * self.t['t3']
        v = acc / 2 * (self.t['t1'] + self.t['t2'] * 2 + self.t['t3'])
        wd += v * self.t['t4']
        wd += v * self.t['t5'] - dec / 6 * (self.t['t5'] ** 2)
        wd += self.t['t6'] * (v - dec / 2 * self.t['t5'] - dec / 2 * self.t['t6'])
        wd += dec / 6 * (self.t['t7'] ** 2)

        return wd


class PatternFrame:
    def __init__(self, src_path):
        # 输入必须是SMEC_motor_a.xlsx的格式
        try:
            df_low_spd = pd.read_excel(src_path, sheetname=1, skiprows=3)
            df_high_spd = pd.read_excel(src_path, sheetname=2, skiprows=3)
            self.df = pd.concat([df_low_spd, df_high_spd])[['MTYPE', 'SPD']]
        except:  # todo: 确认错误类型
            print(u"作为曳引机型号和速度的输入表格的路径或格式有误！")

    def frame(self, method, writer, **kwargs):
        df = self.df.copy().drop_duplicates()
        df.index = range(len(df))
        # print(df)
        try:
            if method == "distance":
                distance = kwargs['distance']
                pat_segs = df.apply(lambda s: Pattern(s.MTYPE, s.SPD).calc_by_distance(distance), axis=1)
                label = 'distance_{0:.2f}'.format(distance)
                # print(pat_segs)
            elif method == 'tr':
                tr = kwargs['tr']
                lr = kwargs['load_ratio']
                df['tr'] = tr
                df['load_ratio'] = lr
                pat_segs = df.apply(lambda s: Pattern(s.MTYPE, s.SPD).calc_by_tr(tr, lr), axis=1)
                label = 'tr_{0:.1f}_{1:d}%'.format(tr, int(lr * 100))
            elif method == 'sph':
                sph = kwargs['sph']
                lr = kwargs['load_ratio']
                df['sph'] = sph
                df['load_ratio'] = lr
                pat_segs = df.apply(lambda s: Pattern(s.MTYPE, s.SPD).calc_by_sph(sph, lr), axis=1)
                label = 'sph_{0:d}_{1:d}%'.format(sph, int(lr * 100))
            else:
                print("速度曲线计算方法中没有" + str(method))
            pat_segs = pd.DataFrame(list(pat_segs.values))  # 若不用list转成list类型，则会保持ndarray类型，没法df
            pat_segs.columns = ['acc', 'dec'] + ['t%i' % i for i in np.arange(1, 9)] + ['wd']
            df = df.join(pat_segs)
        except KeyError:
            print("method对应的变量参数名有误！")

        if writer:
            try:
                df.to_excel(writer, sheet_name=label)
            except:
                pass

        return df


if __name__ == "__main__":
    # src_path = ur"D:\Prjts\2012_CalTables\00.Ref\003.部件整理\SMEC_Motors_a.xlsx"
    src_path = ur"D:\Prjts\2012_CalTables\00.Ref\000.系统计算用表\SMEC_Motors_a.xlsx"
    pf = PatternFrame(src_path)
    writer = pd.ExcelWriter('patterns.xlsx')
    pf.frame(method='distance', writer=writer, distance=3.0)
    pf.frame(method='tr', writer=writer, tr=100, load_ratio=.5)
    pf.frame(method='sph', writer=writer, sph=180, load_ratio=.5)
    writer.save()
