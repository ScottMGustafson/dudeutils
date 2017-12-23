import numpy as np

test_xml = """
    <?xml version="1.0"?>
    <SpecTool version="1.0">
        <CompositeSpectrum id="/home/scott/research/J0744+2059/order_one.fits">
        <Spectrum spec="order_one.fits" error="order_one_e.fits" chi2="3403.764373117464" pixels="1475.0" params="13"/>
        <ContinuumPoint x="3145.012755" y="-3.511155E-16" xError="0.0" yError="0.0" xLocked="true" yLocked="true"
                        id="null"/>
        <ContinuumPoint x="3145.012755" y="-3.511155E-16" xError="0.0" yError="0.0" xLocked="true" yLocked="true"
                        id="test_pt"/>
        <Absorber ionName="C IV" N="12.906365" b="6.010237" z="1.921971" NError="0.0" bError="0.0" zError="0.0"
                  NLocked="false" bLocked="false" zLocked="false" id="CIV7"/>
        </CompositeSpectrum>
        <Region start="3234.068004" end="3235.450664"/>
        <Region start="3235.0" end="3237.322969"/>
        <Region start="3239" end="3239.408987"/>
        <Region start="3247.901669" end="3248.086395"/>
    </SpecTool>
    """.strip()


class FakeModel(object):
    def __init__(self, **kwargs):
        for key, val in kwargs.items():
            setattr(self, key, val)
        self.absorber_list = [MockAb(id='ab1'), MockAb(id='ab2'), MockAb(id='ab3'),
                              MockAb(id='ab4'), MockAb(id='ab5')]
        self.cont_point_list = [MockCont(id='cont1'), MockCont(id='cont2'),
                                MockCont(id='cont3'), MockCont(id='cont4'),
                                MockCont(id='cont5')]
        self.region_list = [MockRegion(0., 10000.)]
        self.spec = MockSpec()

    def set_absorbers(self, ab_lst):
        for ab in ab_lst:
            for i, ab_ in enumerate(self.absorber_list):
                if ab_.id == ab.id:
                    self.absorber_list[i] = ab
                    break

    def set_cont_points(self, cont_lst):
        for cnt in cont_lst:
            for i, cnt_ in enumerate(self.cont_point_list):
                if cnt_.id == cnt.id:
                    self.cont_point_list[i] = cnt
                    break


class MockSpec(object):
    def __init__(self, **kwargs):
        for key, val in kwargs.items():
            setattr(self, key, val)
        self.waves = np.linspace(3000, 5000, num=1000)
        self.flux = np.random.normal(loc=10., scale=1., size=1000)
        self.error = np.sqrt(self.flux)
        self.hdu = None
        self.abs, self.cont = 10. * np.ones(1000), 10. * np.ones(1000)


class MockCont(object):
    def __init__(self, **kwargs):
        self.id = kwargs.get('id', 'contid%d' % np.random.randint(0, 2 ** 32))
        self.x = kwargs.get('x', np.fabs(np.random.normal(4000, 3000)))
        self.y = kwargs.get('y', np.fabs(np.random.normal(10E-14, 10E-15)))
        self.xLocked = kwargs.get('xLocked', False)
        self.yLocked = kwargs.get('yLocked', False)
        for key, val in kwargs.items():
            setattr(self, key, val)


class MockAb(object):
    def __init__(self, **kwargs):
        self.id = kwargs.get('id', 'abid%d' % np.random.randint(0, 2 ** 32))
        self.N = kwargs.get('N', np.fabs(np.random.normal(15.5, 2.)))
        self.b = kwargs.get('b', np.fabs(np.random.normal(16., 2.)))
        self.z = kwargs.get('z', np.fabs(np.random.normal(2.93, 0.05)))
        self.NLocked = kwargs.get('NLocked', False)
        self.bLocked = kwargs.get('bLocked', False)
        self.zLocked = kwargs.get('zLocked', False)
        self.ionName = kwargs.get('ionName', 'H I')

        for key, val in kwargs.items():
            setattr(self, key, val)

    def __eq__(self, other):
        for k in "id N b z ionName".split(" "):
            if getattr(self, k) != getattr(other, k):
                return False
        return True

    def __ne__(self, other):
        return not self == other

    def get_lines(self):
        return [MockLine() for _ in range(np.random.randint(1, 17))]

    def in_region(self, regions):
        """
        detects whether of not CENTER of line is in a optimization region.
        only needs to be in one region to pass
        """
        for region in regions:
            for line in self.get_lines():
                if line.get_obs(z=self.z) in region:
                    return True
        return False


class MockLine(object):
    def __init__(self, **kwargs):
        self.ionName = kwargs.get('ionName', 'H I')
        self.wave = kwargs.get('wave', np.random.uniform(900., 1500.))
        self.f = kwargs.get('f', np.random.uniform(0.01, 1.5))
        self.rest_wave = self.wave

    def get_obs(self, z):
        return (1. + z) * self.wave


class MockRegion(object):
    def __init__(self, start, end):
        self.start, self.end = start, end

    def __contains__(self, item):
        return self.start <= item <= self.end

    def __gt__(self, other):
        return self.end > other.end

    def __lt__(self, other):
        return self.start < other.start

    def __ge__(self, other):
        return self.end >= other.end

    def __le__(self, other):
        return self.start <= other.start
