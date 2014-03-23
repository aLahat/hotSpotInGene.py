import Image
import pickle
import os
try:import winsound
except:pass

class chrImg():
# read and make image
    w = 1000 
    def __init__(self,CHR):
        try: self.chrImg = Image.open(CHR + '.jpg')
        except : 
            self.CHR = CHR
            try:
                d = pickle.load( open( "save.p", "rb" ) ) 
                self.d = d[CHR]
            except:
                print 'reading GTF'
                d = self.makeGTFdict('Mus_musculus.GRCm38.74.gtf')
                pickle.dump(d, open( "save.p", "wb" ) )
                self.d = d[CHR]
            length = 0
            for gene in self.d:
                if length < self.d[gene][1]: length=self.d[gene][1]
            self.length = length
            self.h = length/self.w + 1
            self.chrImg = Image.new('RGB', (self.w,self.h), (0,0,0))
            self.colorChr()
 
    def addBase(self, base, strand):
        coords = self.baseCoords(base)
        currCol = self.chrImg.getpixel(coords)
        if strand == '+': col = (currCol[0]+1, currCol[1], currCol[2])
        if strand == '-': col = (currCol[0], currCol[1]+1, currCol[2])
        self.chrImg.putpixel(coords,col)
        
    def colorChr(self):
        for gene in self.d:
            start= self.d[gene][0]
            end= self.d[gene][1]
            strand= self.d[gene][2]
            for base in range(start,end):
                self.addBase(base,strand)
#        self.chrImg.save(self.CHR + '.jpg')
        
    def getGTFline(self, line): 
        details = line.split('\t')
        name=details[-1].split('gene_name "')[-1]
        name = name[:name.find('"')]
        chromosome =details[0]
        TYPE =details[1]
        dictionary={'chr':chromosome,'type':TYPE,'start':details[3],'end':details[4],'strand':details[6],'name':name}
        return dictionary   
    def makeGTFdict(self, FILE, allowedTypes = ['protein_coding']):
        f= open(FILE,'r')
        bothDicts={}
        names=[]
        while True:
            line = f.readline()
            if line == '':
                break
            line = self.getGTFline(line)
            if not(line['name'] in names):
                names.append(line['name'])
                if not line['chr'] in bothDicts: 
                    bothDicts.update({line['chr']:{}})
                    print line['chr']
                if line['name'] in bothDicts[line['chr']]:
                    if int(line['start']) < bothDicts[line['chr']][line['name']][0]:
                        bothDicts[line['chr']][line['name']][0]=int(line['start'])
                    if int(line['end']) > bothDicts[line['chr']][line['name']][1]:
                        bothDicts[line['chr']][line['name']][1]=int(line['end'])
                else:
                    bothDicts[line['chr']].update({line['name']:[int(line['start']),int(line['end']),line['strand']]})
            
        f.close()
        return bothDicts        
    def baseCoords(self, base):
        w = self.w
        x = base%w
        y = base/w
        return (x,y)

    def ingeneBin(self, steps,binSize):
        output = dict.fromkeys(range(0, self.length, steps),0)
        stepD = dict.fromkeys(range(0,self.length,steps))
        n=0
        for step in stepD:
            if n%100==0: print 'making stepD: '+str(int(float(n)/len(stepD)*100))+'%'
            start = step
            end = start + steps
            value = 0
            prev=0
            for base in range(start,end):
                try:
                    coords = self.baseCoords(base)
                    genes = sum(self.chrImg.getpixel(coords))
                    if genes>1 and not genes == prev:
                        value+=1
                        prev = genes
                except: pass 
                stepD[step] = value  
            n+=1
        self.stepD = stepD
        n=0 
        
        #############################################################
        for step in output.keys():
            if n%100==0: print 'making output: '+str(int(float(n)/len(output)*100))+'%'
            end =  binSize+ step
            value = 0
            for i in range(step, end, steps):
                try: output[step] += stepD[i]
                except: pass
            n+=1
        self.output = output
        txt = []
        for i in output:
            line = str(i)+'\t'+str(output[i])+'\n'
            txt.append(line)
        txt=''.join(txt)
        f = open(self.CHR+'.geneHistogram.csv','w')
        f.write(txt)
        f.close()
    def strandBin(self, steps,binSize):
        output = dict.fromkeys(range(0, self.length, steps),0)
        for i in output:     
            pass


print 'start'
allCHR = ['JH584293.1',
 'GL456354.1',
 'MG4180_PATCH',
 'GL456219.1',
 'GL456381.1',
 'MG132_PATCH',
 'MG4136_PATCH',
 'GL456385.1',
 'JH584298.1',
 'MG4213_PATCH',
 'MG153_PATCH',
 'GL456372.1',
 'MG4151_PATCH',
 '3',
 '2',
 '5',
 '4',
 '7',
 'MG4211_PATCH',
 '9',
 'MG4222_MG3908_PATCH',
 'GL456211.1',
 'JH584296.1',
 'MG4209_PATCH',
 'MG3829_PATCH',
 'GL456221.1',
 'GL456350.1',
 '6',
 'GL456233.1',
 'MG4212_PATCH',
 '17',
 'MG4214_PATCH',
 '8',
 'Y',
 'X',
 'GL456239.1',
 'JH584294.1',
 '11',
 '10',
 '13',
 '12',
 '15',
 '14',
 'GL456210.1',
 '16',
 '19',
 '18',
 'JH584292.1',
 'JH584295.1',
 'GL456216.1',
 'MT',
 'JH584297.1',
 '1',
 'GL456212.1',
 'JH584304.1',
 'MG3833_PATCH',
 'JH584299.1',
 'JH584303.1']
for CHR in allCHR:
    print '\n'+CHR+'\n'
    c = chrImg(CHR)
    c.ingeneBin(100000,500000)
try: winsound.Beep(1000,1000)
except:pass
#im = Image.new('RGB', size, 'white')
#im.save('1000_1000.jpg')
