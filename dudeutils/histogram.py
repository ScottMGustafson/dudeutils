import matplotlib.pyplot as plt
import numpy as np


plt.rc('text', usetex=True)
plt.rc('font', family='serif', size=18)


def histogram(data, num_bins=20, density=False):
    """plots a histogram.  If you need any filtering of data, do that before hand"""

    hist, _x = np.histogram(data, bins=num_bins, density=False)

    x=[]    
    for i in range(len(hist)):
        x.append((_x[i]+_x[i+1])/2.)

    if density:
        hist = normalize_hist(hist)

    return x, hist


def plot_histogram(data, num_bins=20, xlabel=None, ylabel=None, title=None, density=False):
    """plots a histogram.  If you need any filtering of data, do that before hand"""

    hist, _x = np.histogram(data, bins, density)

    x=[]    
    for i in range(len(hist)):
        x.append((_x[i]+_x[i+1])/2.)

    if density:
        hist = normalize_hist(hist)

    x, hist = histogram(data, num_bins, xlabel=None, ylabel=None, title=None, density=False)

    fig,ax = plt.subplots(figsize=(9.,7.))

    plt.plot(x, hist,linestyle='steps')
    try:
        plt.title(title.replace('_','-'))
    except:
        plt.title('')
    ax.get_xaxis().get_major_formatter().set_useOffset(False)
    ax.get_yaxis().get_major_formatter().set_useOffset(False)
    ax.set_xlabel(xlabel)
    ax.set_ylabel(ylabel)
    ax.set_ylim([0.,(1.05)*max(hist)])
    plt.subplots_adjust(top=0.87,bottom=0.13,left=0.1,right=0.95)
    #plt.savefig(output_fname)    
    print('n=',len(data))
    plt.show()

def normalize_hist(hist):
    tot = 0
    for item in hist:
        tot+=float(item)
    return [item/tot for item in hist]

def bootstrap(data):
    data=np.array(data)
    indices = np.random.randint(len(data),size=len(data))
    return data[indices]

def bootstrap_resampling(data,samples=1000, num_bins=20, xlabel=None, ylabel=None, title=None, density=False):
    #total = []
    bootstrap_dist = []
    for i in range(samples):
        shuffled = bootstrap(data)
        bootstrap_dist.append(np.mean(shuffled))
        #total=np.concatenate( (total, shuffled, axis=0 )
    print('n=',len(data),'SE = ',np.std(bootstrap_dist),'mu = ',np.mean(bootstrap_dist),'med = ',np.median(bootstrap_dist),1.96*np.std(bootstrap_dist)*np.sqrt(len(data)))
    histogram(bootstrap_dist, num_bins, xlabel, ylabel, title, density)
    return bootstrap_dist

if __name__ == '__main__':
    data=np.array([-4.5241548084476335, -4.525209059705556, -4.529532905915431, -4.532590409137978, -4.532586427400982, -4.519782990779106, -4.519782103503296, -4.5352731636484265, -4.550804481835167, -4.551395321855489, -4.582699666398964, -4.646214672128641, -4.76639688005803, -4.738284515111786, -4.595131253375801, -4.594162780113242, -4.780837575794598, -4.622991714281932, -4.689698554847919, -4.64493449269642, -4.511083300767959, -4.483604507471817, -4.533189294682931, -4.501966908493383, -4.499224675851444, -4.51706572263959, -4.504183579373066, -4.646830854495326, -4.569774663075371, -4.5486296391886984, -4.521851613510405, -4.529968206515765, -4.48538096375643, -4.485348495133719, -4.493392607380926, -4.491670177781511, -4.4893369127272535, -4.489995350799225, -4.489546528002478, -4.487301324275409, -4.487536580875068, -4.487892296104173, -4.523627446980871, -4.521903990498432, -4.516347015769043, -4.529378227511888, -4.5481354441795645, -4.547665200256468, -4.560839710622224, -4.505571530692524, -4.5064137568969205, -4.508090319187783, -4.509043628383498, -4.499977760478458, -4.573392199999999, -4.570302331908625, -4.580809191567784, -4.554917531363891, -4.554602829957686, -4.548264988595699, -4.548753563140254, -4.549596612914106, -4.549253894846448, -4.551938423653205, -4.550180743990083, -4.543105044564422, -4.532724485267371, -4.542599199381874, -4.551236602638804, -4.551815777825313, -4.56262375506066, -4.548672996240592, -4.5479221272450285, -4.547894165595759, -4.571942774165864, -4.595476009238178, -4.527295465580426, -4.557537283124834, -4.557576920237846, -4.541055766044623, -4.5535050865904925, -4.569122132068381, -4.5643764462418055, -4.563253288144962, -4.562139435523703, -4.525721715443412, -4.525721028692498, -4.532807017992663, -4.732458495825414, -4.584498958499006, -4.583573313457972, -4.5630494416623755, -4.560996974143199, -4.560592037715693, -4.557157460319948, -4.556953065135517, -4.555198030792667, -4.535715542306001, -4.535714932377145, -4.560006477221249, -4.576888488759822, -4.593503867891686, -4.5670384246043145, -4.554934455743313, -4.542893844829804, -4.542882150233444, -4.530770816455915, -4.555527817412299, -4.554942624619368, -4.560384768895336, -4.565472291942363, -4.610068315870784, -4.575714966874326, -4.561531529561297, -4.5615368291069185, -4.570802617924535, -4.567683505215282, -4.566387720810125, -4.581499874937428, -4.548614903417192, -4.621260048903224, -4.609316425867295, -4.6634897960847645, -4.605835119009864, -4.584391629007845, -4.565589838004174, -4.600711741046403, -4.598870050090218, -4.565482754188247, -4.565482754188247, -4.545941985353659, -4.571854627230806, -4.569353623916342, -4.5873768998328845, -4.63729997093165, -4.636705140109624, -4.636139496759698, -4.654004297550774, -4.503346865853404, -4.550475108753197, -4.495964981716654, -4.496332597686212, -4.509389770610467, -4.508814131612054, -4.510738452945093, -4.510807432422112, -4.5102247071203045, -4.501976015004772, -4.502447199942292, -4.503209411597982, -4.518236085504613, -4.516606745220466, -4.510936048105842, -4.515052640857014, -4.512404355453853, -4.5122535506485, -4.5122535506485, -4.5565484540880945, -4.533409607887316, -4.533404271963066, -4.5151928986975935, -4.515186279947491, -4.552798215430055, -4.515377108574539, -4.552872408803012, -4.506203630810086, -4.529450991179228, -4.552576279399135, -4.553102867924485, -4.529401851318237, -4.539766010612725, -4.937586600308821, -4.701381324959696])

    bootstrap_resampling(data, samples=6700)

