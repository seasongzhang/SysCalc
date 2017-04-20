# -*- coding:utf-8 -*-

from __future__ import division


class Pattern:
    def __init__(self, mtype, spd, method):
        """

        :param mtype: 曳引机型号，如"ZPML-J117","ZPML-GH70"
        :param spd: 额定速度，单位为m/min
        """

        self.spd = float(spd)
        self.method = method

        power = mtype.split('-')[1][1:4]
        try:
            is_low_spd = int(power) < 400
        except ValueError:
            is_low_spd = False

        if not is_low_spd:
            self.acc, self.dec = [0.8, 0.8]
            self.t = {'t1': 1.1, 't3': 1.0, 't5': 1.0, 't7': 1.2}
        elif self.spd >= 150:
            self.acc, self.dec = [0.8, 0.8]
            self.t = {'t1': 0.8, 't3': 0.8, 't5': 0.8, 't7': 0.8}
        else:
            self.acc, self.dec = [0.6, 0.6]
            self.t = {'t1': 0.6, 't3': 0.6, 't5': 0.6, 't7': 0.6}

    def calc_by_tr(self, tr=150, load_ratio=0.5):
        """

        :param tr: 提升高度
        :param load_ratio: 负载持续率，根据提升高度计算得到电梯单次运行时间后根据负载持续率计算休息时间
        :return:
        """
        v = self.spd / 60
        self.t['t2'] = v / self.acc - (self.t['t1'] + self.t['t3']) / 2
        self.t['t6'] = v / self.dec - (self.t['t5'] + self.t['t7']) / 2
        s = [self.acc / 6 * self.t['t1'] ** 2,  # t1内的行程
             self.acc / 2 * (self.t['t1'] + self.t['t2']) * self.t['t2'],  # t2内的行程
             self.acc * (self.t['t1'] / 2 + self.t['t2'] + self.t['t3'] / 3) * self.t['t3'],  # t3内的行程
             v * self.t['t5'] - self.dec / 6 * self.t['t5'] ** 2,  # t5内的行程
             self.t['t6'] * (v - self.dec / 2 * self.t['t5'] - self.dec / 2 * self.t['t6']),  # t6内的行程
             self.dec / 6 * self.t['t7'] ** 2  # t7内的行程
             ]

        # 必须判断tr>s
        s_uniform = tr - sum(s)
        self.t['t4'] = s_uniform / v
        t_total = sum(self.t.values()) / load_ratio
        self.t['t8'] = t_total * (1 - load_ratio)
        # print(sum(s))

    def calc_by_distance(self, distance):
        """
        
        :param distance: 
        :return: 
        根据在Mathmatica中的计算，针对唐启峰提供的3套加速度/减速度/圆角时间的参数，若仅保证4段圆角时间，3套参数对应的
        行程分别为0.432m/1.024m/1.736m。因此，只要保证层间距大于1.736m，单次行程就必然包含完整的4个圆角，这在一定程
        度上简化了模型（计算）。
        """
        if distance < 1.736:  # todo: 针对不同的梯种设置不同的限制
            print(u'不允许单次行程小于1.736m！')
        else:
            v = self.spd / 60
            self.t['t2'] = v / self.acc - (self.t['t1'] + self.t['t3']) / 2
            self.t['t6'] = v / self.dec - (self.t['t5'] + self.t['t7']) / 2
            s = [self.acc / 6 * self.t['t1'] ** 2,  # t1内的行程
                 self.acc / 2 * (self.t['t1'] + self.t['t2']) * self.t['t2'],  # t2内的行程
                 self.acc * (self.t['t1'] / 2 + self.t['t2'] + self.t['t3'] / 3) * self.t['t3'],  # t3内的行程
                 v * self.t['t5'] - self.dec / 6 * self.t['t5'] ** 2,  # t5内的行程
                 self.t['t6'] * (v - self.dec / 2 * self.t['t5'] - self.dec / 2 * self.t['t6']),  # t6内的行程
                 self.dec / 6 * self.t['t7'] ** 2  # t7内的行程
                 ]
            s_uniform = distance - sum(s)
            if s_uniform >= 0:
                s.insert(index=3, object=s_uniform)
                self.t['t4'] = s_uniform / v
            else:
                pass

    def calc_by_sph_3phase(self, sph=180, load_ratio=0.5):
        """
        由sph和load_ratio可以决定工作制
        :param sph: start per hour
        :param load_ratio: 负载持续率，无论sph值如何，保证每小时内有load_ratio%的时间电梯处于非停止状态
        :return:
        """
        self.t['t2'] = self.spd / 60 / self.acc - (self.t['t1'] + self.t['t3']) / 2
        self.t['t6'] = self.spd / 60 / self.dec - (self.t['t5'] + self.t['t7']) / 2

        t4 = load_ratio * 3600 / sph - sum(self.t.values())

        # 若计算得到的匀速时间>0，说明按照sph和负载持续率计算确实存在匀速运行阶段
        # 若计算得到的匀速时间<0，说明按照sph和负载持续率计算不存在匀速运行阶段，需要考虑加速过程中直接减速的情况
        if t4 >= 0:
            self.t['t4'] = t4
        else:
            self.t['t4'] = 0
            t2_plus_t6 = 3600 / sph - (self.t['t1'] + self.t['t3'] + self.t['t5'] + self.t['t7'])
            t2_minus_t6 = ((self.t['t5'] + self.t['t7']) - (self.t['t1'] + self.t['t3'])) / 2
            self.t['t2'] = (t2_plus_t6 + t2_minus_t6) / 2
            self.t['t6'] = (t2_plus_t6 - t2_minus_t6) / 2
        self.t['t8'] = 3600 / sph * (1 - load_ratio)

    def clac_by_method(self):
        if self.method == 'sph_3phase':
            self.calc_by_sph_3phase()
        elif self.method == 'tr':
            self.calc_by_tr()


if __name__ == "__main__":
    pat = Pattern(mtype='ZPML-GH70', spd='480', method='sph_3phase')
    pat.calc_by_sph_3phase(sph=180, load_ratio=0.5)
    for k in sorted(pat.t.keys()):
        print(k, pat.t[k])
    pat.calc_by_tr(tr=250, load_ratio=0.5)
    for k in sorted(pat.t.keys()):
        print(k, pat.t[k])
