import numpy as np
import matplotlib.pyplot as plt
from scipy.interpolate import griddata
class ParseSplatData(object):
    '''

    '''
    def __init__(self, ref_lat, ref_lon, limit_x, limit_y, grid_x, grid_y ):
        '''

        '''
        self.ref_lat = ref_lat
        self.ref_lon = ref_lon
        self.limit_x = limit_x
        self.limit_y = limit_y
        self.grid_x = grid_x
        self.grid_y = grid_y

    # def haversine_distance(self, ref_loc, locs):
    #     '''
    #
    #     :param reference:
    #     :param locs:
    #     :param earth_radius_km:
    #     :return:
    #     '''
    #     #using Haversine formula
    #     earth_radius_km = 6371
    #
    #     diff = np.radians( locs - ref_loc )
    #
    #     dlat, dlon = diff[:, 0], diff[:, 1]
    #     ref_loc_lat, locs_lat = ref_loc[:, 0], locs[:, 0]
    #     a = np.power( np.sin(dlat / 2), 2.0) + np.cos( np.radians(ref_loc_lat)) * np.cos(np.radians(locs_lat)) * np.power( np.sin(dlon / 2) , 2.0)
    #     c = 2 * np.arcsin( np.sqrt(a))
    #     d = c * earth_radius_km *1000
    #
    #     return d

    def getXCoords(self, ref_lon, locs):
        '''

        :param ref_lon:
        :param locs:
        :return:
        '''
        earth_radius_km = 6371

        locs_lat = locs[:, 0]
        locs_lon = locs[:, 1]

        locs_lat_rad = np.radians( locs_lat )
        dlon_rad = np.radians(locs_lon - ref_lon)

        a =  np.power( np.cos( locs_lat_rad ), 2.0)  * np.power( np.sin(dlon_rad / 2) , 2.0)

        c = 2 * np.arcsin( np.sqrt(a))
        d = c * earth_radius_km *1000

        polarized_d = np.where( locs_lon < ref_lon, -d, d )

        return polarized_d

    def getYCoords(self, ref_lat, locs):
        '''

        :param ref_lat:
        :param locs:
        :return:
        '''
        earth_radius_km = 6371

        locs_lat = locs[:, 0]
        #locs_lon = locs[:, 1]

        dlat_rad = np.radians( locs_lat - ref_lat )

        a = np.power(np.sin(dlat_rad / 2), 2.0)
        c = 2 * np.arcsin( np.sqrt(a))
        d = c * earth_radius_km *1000

        polarized_d = np.where( locs_lat < ref_lat, -d, d )

        return polarized_d

    def convertToEuclidian(self, locs):
        '''

        :param ref_lat:
        :param ref_lon:
        :param locs:
        :return:
        '''
        xcoords = self.getXCoords(self.ref_lon, locs)
        #print "min x-coord, max x-coord: ", np.min(xcoords), np.max(xcoords)
        ycoords = self.getYCoords(self.ref_lat, locs)
        #print "min y-coord, max y-coord: ", np.min(ycoords), np.max(ycoords)
        coords = list(zip(xcoords, ycoords))
        return coords


    def generateInterpolatedMap(self, coords, vals, outfile, showMap = False):
        '''

        :param coords:
        :param vals:
        :return:
        '''

        xi = np.linspace( -self.limit_x, self.limit_x, self.grid_x )
        yi = np.linspace( -self.limit_y, self.limit_y, self.grid_y )
        grid_x, grid_y =  np.meshgrid( xi, yi )
        #zi = griddata(coords[:, 0], coords[:, 1], vals, xi, yi, interp='linear')
        zi_with_nan = griddata(coords, vals, (grid_x, grid_y), method='cubic')
        zi_nearest = griddata(coords, vals, (grid_x, grid_y), method='nearest')
        zi = np.where(np.isnan(zi_with_nan), zi_nearest, zi_with_nan)
        np.savetxt(outfile+".txt", zi, fmt= '%.4f' )

        if showMap:
            #CS = plt.contour(xi, yi, zi, 10, linewidths= 0.5, colors='k')
            plt.contourf(xi, yi, zi, 10, cmap = 'viridis')
            plt.colorbar()
            plt.savefig(outfile+".png")
            plt.close()



    def filterCoords(self, coords, vals):
        '''

        :return:
        '''
        scoords, svals = [], []

        for c, v in zip(coords, vals):
            if (-self.limit_x <= c[0] <= self.limit_x) and (-self.limit_y <= c[1] <= self.limit_y):
                scoords.append( (c[0], c[1] ) )
                svals.append( v )
        return np.array(scoords), np.array(svals)


    def generateMap(self, inputSplatFileName, outputMapFileName ):
        '''

        :param fileName:
        :return:
        '''
        lat_lon_list, val_list = [], []
        with open(inputSplatFileName, 'r') as f:
            _, _ = f.readline(), f.readline()
            for line in f:
                line_vals =  line.split(',')
                lat , lon, val = float(line_vals[0]), -float(line_vals[1]), float( line_vals[4].split()[0] )
                lat_lon_list.append((lat, lon))
                val_list.append(val)

            vals = np.array(val_list, dtype=np.float)
            locs = np.array(lat_lon_list, dtype= np.float)

            coords = self.convertToEuclidian(locs)
            scoords, svals = self.filterCoords(coords, vals)
            self.generateInterpolatedMap(scoords, svals, outputMapFileName, showMap=True)



if __name__ == '__main__':
    np.random.seed(1009993)
    print ("Parsing Splat Data files")
    splatFileName = 'pathloss_2.dat'
    mapFileName = 'pathloss_2'
    ref_lat, ref_lon = 40.890308, -73.127174
    limit_x, limit_y = 6000.0, 6000.0

    grid_x, grid_y = 100, 100
    pd = ParseSplatData(ref_lat, ref_lon, limit_x, limit_y, grid_x, grid_y)
    pd.generateMap(splatFileName, mapFileName)
