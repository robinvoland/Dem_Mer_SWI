import numpy as np
import flopy as fp
import geopandas as gpd
import pickle
import rasterstats
import rasterio
import pandas as pd
import scipy as sp
import os


def gp2cellids (grid, gdf, type = "polygon",layer=0,a_fac=3):
    
    """
    this function extract the cellids of the intersection between a geopandas object and a grid 
    grid : modelgrid with flopy.discretisation !
    gdf : geopandas object with one entity only
    idomain : array, the idomain array to update it
    idomain_active : bool, if true the idomain is update (cells intersect by the gdf will be noted as active), prevents some issues
    type : str, features type (polygon or line)
    layer : int, the layer on which is the gdf
    areas : factor that determine if a cell is accounted as intersected or not based on the total area intersected
    (a value of 3, for example, means only cells which have 1/3 of their area intersected by the polygon will be taken into account)
    """
    
    ix = fp.utils.gridintersect.GridIntersect(grid)
    if type == "polygon":
        result = ix.intersect(gdf.geometry[0])
        result = result[result.areas>(np.nanmax(result.areas)/a_fac)] # only take into account cells that have a least 1/3 intersected 
        
        
    if type == "boundary" :
        result = ix.intersect(gdf.geometry[0].boundary)
        
    if type == "line" :
        result = ix.intersect(gdf.geometry[0])
        
    result = result[result.areas!=0]                       # fix bug with some null areas
    
    lst=[]
    for irow, icol in result.cellids:
        lst.append(((layer,irow,icol)))

    return lst

def write_pickle(path, data):
    '''
    Function to write a picke file.
    '''
    
    with open(path,'wb')as f:
        pickle.dump(data, f, pickle.HIGHEST_PROTOCOL)
        
def read_pickle(path):
    '''
    Function to read a pickle file.
    '''
    
    with open(path, 'rb')as file:
        f_read = pickle.load(file)
    
    return f_read

def rspl_rast(rast_path,grid,band=1):
    
    """
    Use the resample_to_grid method from flopy Raster. 
    rast_path : path to the raster
    grid : modelgrid (gwf.modelgrid or flopy.discretisation)
    """
    
    rast = fp.utils.Raster.load(rast_path)
    arr = rast.resample_to_grid(grid,band)
    return arr


def liss_mob(arr,n,null_v = 0):
    
    """
    Apply a moving average (with 2*n numbers) on 2D array.
    arr : 2D numpy array
    n : number of elements (in one of the four direction) to take into account for the moving average (n=2 --> average of a specific number will be calculated with the surroundings 5x5 elements)
    return a 2D array and replace null value by 0
    """
    
    
    arr[arr==null_v]=None
    for irow in range(n,arr.shape[0]-n):
        for icol in range(n,arr.shape[1]-n):
            if not np.isnan(arr[irow,icol]):
                bloc = arr[irow-n:irow+n+1,icol-n:icol+n+1]
                arr[irow,icol] = np.nanmean(bloc)
    arr = np.nan_to_num(arr)
    return arr


def import_riv(Riv_path,MNT_path,data_dir='./',max_length=20):
    
    MNT = rasterio.open(MNT_path)
    Riv = gpd.read_file(Riv_path)
    
    Riv_seg = Riv.geometry.segmentize(max_segment_length=max_length) #Cut the rivers into lines of small lenght
    points=gpd.GeoSeries.from_xy(Riv_seg.get_coordinates().x,
                                 Riv_seg.get_coordinates().y) #retrieve the coordinates of every vertices
    alts = rasterstats.point_query(points.unary_union, 
                                      MNT_path,
                                      interpolate='nearest') #get the value of the MNT at those points
    gdf=gpd.GeoDataFrame(geometry=points).reset_index(drop=True) #create a geodataframe with all this points
    gdf['alt']=alts[0]
    
    gdf.to_file(data_dir+Riv_path.split('/')[-1][:-9]+'.shp')
    
    
def get_Total_Budget(model_name,model_dir):

    """
    Return a DF containing Budget data for the entire model in the LST file
    Only for the 1st time step, 1st stress period for the moment
    model_name : str, name of the model given in the gwf pack
    model_dir : str, path to workspace
    npack : number of additionnal packages with save_flows set True
    """
    
    
    file = "{}/{}.lst".format(model_dir,model_name)
    f = open(file,"r")
    i=-1
    for ilin in f.readlines():
        i += 1
        if ilin =='  VOLUME BUDGET FOR ENTIRE MODEL AT END OF TIME STEP    1, STRESS PERIOD   1\n': # check at which line the budget is
            break
    
    ###number of packages
    npack=0
    for o in range(100):
        f = open("{}/{}.lst".format(model_dir,model_name),"r")
        if f.readlines()[i+8+o]=="\n":
            break
        npack +=1
    ###number of packages
    
    # retrieve data
    lst_val_IN =[]
    lst_val_OUT = []
    lst_nam_pak = []
    pak_type=[]
    for ipak in range(npack):
        ipak += 8
        
        f = open("{}/{}.lst".format(model_dir,model_name),"r")
        lst_nam_pak.append(f.readlines()[i+ipak][85:96].rstrip())

        f = open("{}/{}.lst".format(model_dir,model_name),"r")
        lst_val_IN.append(float(f.readlines()[i+ipak][63:80]))

        f = open("{}/{}.lst".format(model_dir,model_name),"r")
        lst_val_OUT.append(float(f.readlines()[i+ipak+npack+5][63:80]))
        
        f = open("{}/{}.lst".format(model_dir,model_name),"r")
        pak_type.append(f.readlines()[i+ipak][58:62])

    Budget = pd.DataFrame({"Pack":lst_nam_pak,
                  "IN":lst_val_IN,
                 "OUT":lst_val_OUT,
                  "Type":pak_type})

    return Budget


def budget_vert(gwf,head,model_name,model_dir,Liste_epontes,name_epontes):
    cbb = fp.utils.CellBudgetFile(os.path.join(model_dir, model_name+".cbc"))
    spdis = cbb.get_data(text="SPDIS")[-1]
    qx, qy, qz = fp.utils.postprocessing.get_specific_discharge(
        spdis, gwf, head=head
    )
    
    res={}
    names=[]
    res_in=[]
    res_out=[]
    
    for i in range(len(Liste_epontes)):
        names.append(name_epontes[i])
        res_in.append(-1*np.sum(qz[i][qz[i]<0])*200*200)
        res_out.append(np.sum(qz[i][qz[i]>0])*200*200)
    res['Name']=names
    res['IN']=res_in
    res['OUT']=res_out
    
    return pd.DataFrame(res)
    
    
def get_Riv_Budget(model_name,model_dir,riv_name):

    file = "{}/{}.lst".format(model_dir,model_name)
    f = open(file,"r")
    i=-1
    for ilin in f.readlines():
        i += 1
        if ilin ==" RIV PACKAGE ({}) FLOW RATES   PERIOD      1   STEP        1\n".format(riv_name): # check at which line the budget is
            break

    cell=[]
    flow=[]
    for o in range(1000):
        f = open(file,"r")
        line=f.readlines()[i+4+o]
        if line==' -----------------------------------------------\n':
            break
        line_split=line.split()
        cell.append(line_split[1][1:-1])
        flow.append(line_split[2])

    return cell,flow

def collect_all_rivs_budget(model_name,model_dir,list_rivs,idomain):
    #Positif = IN Unité en m^3/s
    
    res={}
    l=[]
    x=[]
    y=[]
    Q=[]
    river_name=[]
    for name in list_rivs:
    
        cell,flow=get_Riv_Budget(model_name,model_dir,name)


        for i in range(len(cell)):
            cellid=cell[i].split(',')
            l.append(int(cellid[0]))
            x.append(int(cellid[1]))
            y.append(int(cellid[2]))
            Q.append(float(flow[i]))
            river_name.append(name)

        
    res['Layer']=l
    res['X']=x
    res['Y']=y
    res['Flow']=Q
    res['river name']=river_name

    df=pd.DataFrame(res)
    array=np.zeros(np.shape(idomain[0]))

    for i in range(len(df)):
        array[df.X[i],df.Y[i]]=-df.Flow[i]

    return array



def create_botm(surfs):
    surfs=surfs[::-1]
    bot=surfs[0].copy()
    for i in range(len(surfs)-1):
        mask=np.array(surfs[0]==9999)
        for j in range(i):
            mask*=np.array(surfs[j+1]==9999)
            
        bot[mask]=surfs[i+1][mask]
    mask=(surfs[-1]-bot)<10
    bot[mask]=surfs[-1][mask]-10
    bot[bot==9989]=9999
    botm_plot=bot.copy()
    botm_plot[botm_plot==9999]=np.nan
    
    surfs=surfs[:-1]
    surfs=surfs[::-1]
    surfs.append(bot)
    return surfs,botm_plot
    
    
def remove_cell_iso_2D(idomain,min_size=10):
    res=idomain.copy()
    labels=sp.ndimage.label(idomain)
    label_count=np.zeros(labels[-1]+1)
    for j in range(labels[-1]+1):
        label_count[j]=np.sum(labels[0]==j)
    res[np.isin(labels[0],np.arange(labels[-1]+1)[label_count<min_size]) & np.equal(idomain,1)]=0 #The groups smaller than min_size are set to inactive
    return res     
    
def create_idomain(surfs,thick_min=1,min_size=100):
    nlay,ncol,nrow=np.shape(surfs)-np.array([1,0,0])
    idomain=np.ones((nlay,ncol,nrow))
    
    for i in range(nlay):
        idomain[i][surfs[i]-surfs[i+1]<thick_min]=0
        idomain[i][surfs[i+1]==9999]=0
        idomain[i]=remove_cell_iso_2D(idomain[i],min_size)
        idomain[i][surfs[i]-surfs[i+1]<thick_min]=-1
    return idomain


def create_eponte_2_layers(layer,top_layer,bot_layer,thickness=2,eps=0.1):
    
    res=layer.copy()

    res[layer-bot_layer>thickness/2+eps]-=thickness/2
    layer[top_layer-layer>thickness/2+eps]+=thickness/2
    res_plot=res.copy()
    res_plot[res_plot==9999]=np.nan
    return layer,res,res_plot

def create_eponte_sea(top,sea_bol,bot_layer,thickness=2):
    res=top.copy()
    res[sea_bol==1]-=thickness
    bot_layer[bot_layer>res]=res[bot_layer>res]
    res_plot=res.copy()
    res_plot[res_plot==9999]=np.nan
    
    return bot_layer,res,res_plot
        