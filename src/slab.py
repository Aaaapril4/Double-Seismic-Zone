from scipy.spatial import KDTree
from scipy.interpolate import griddata
from pyproj import Proj
import numpy as np

class Slab:
    def __init__(self, grid, region) -> None:
        self.grid = grid  #[longitude, latitude, depth]
        self.region = region #[min_longitude, max_longitude, min_latitude, max_latitude]
        
    
    def transform(self) -> None:
        center = [(self.region[1]+self.region[0]) // 2, 
                  (self.region[3]+self.region[2]) // 2]
        self.proj = Proj(f"+proj=sterea +lon_0={center[0]} +lat_0={center[1]} +units=km")
        lt = self.proj(longitude=self.region[0], latitude=self.region[3])
        lb = self.proj(longitude=self.region[0], latitude=self.region[2])
        rt = self.proj(longitude=self.region[1], latitude=self.region[3])
        rb = self.proj(longitude=self.region[1], latitude=self.region[2])
        self.region = [min(lt[0], lb[0]), max(rt[0], rb[0]), min(lb[1], rb[1]), max(lt[1], rt[1])]
        for i, g in enumerate(self.grid):
            self.grid[i][0], self.grid[i][1] = self.proj(longitude=g[0], latitude=g[1])


    def interp(self, spacing: float) -> None:
        '''
        interpolate the slab
        @spacing: float
        '''
        grid_x, grid_y = np.mgrid[
            self.region[0]:self.region[1]:spacing, 
            self.region[2]:self.region[3]:spacing
        ]

        grid_z = griddata(self.grid[:, 0:2], self.grid[:, 2], (grid_x, grid_y), method='cubic')

        grid_x = grid_x.ravel()
        grid_y = grid_y.ravel()
        grid_z = grid_z.ravel()

        stacked_grid = np.column_stack((grid_x, grid_y, grid_z))

        mask = ~np.isnan(stacked_grid[:, 2])
        stacked_grid = stacked_grid[mask]
    
    
    def init_tree(self):
        self.tree = KDTree(self.grid)
    

    def get_distance(self, event: list) -> float:
        '''
        get the normalized distance between the event and the slab
        @event: [longitude, latitude, depth]
        @return: distance (positive if deeper, negative if shallower)
        '''
        if hasattr(self, 'proj'):
            event = list(self.proj(event[0], event[1])) + event[2:]
        dist, index = self.tree.query(event)
        projection = self.grid[index]
        if projection[2] < event[2]:
            return dist
        else:
            return -dist


if __name__ == '__main__':
    grid = np.loadtxt('/mnt/home/jieyaqi/code/AlaskaEQ/data/slab_Fan.xyz', delimiter=',')
    s = Slab(grid, [-166, -148, 50, 60])
    s.transform()
    s.interp(1)
    s.init_tree()
    print(s.get_distance([-154, 56, 20]))