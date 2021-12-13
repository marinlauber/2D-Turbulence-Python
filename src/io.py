from netCDF4 import Dataset
import numpy as np

class netCDFwriter(object):

    def __init__(self, name, flow) -> None:

        # the data set that we use to store data
        self.data = Dataset(name+'.nc','w','NETCDF4') # using netCDF4 for output format
        
        # two space dimension and a time dimension
        self.data.createDimension('x',flow.nx)
        self.data.createDimension('y',flow.ny)
        self.data.createDimension('t',None)

        # fill-in the coordinates
        self.x = self.data.createVariable('x','float32',('x'))
        self.x[:] = flow.x
        self.y = self.data.createVariable('y','float32',('y'))
        self.y[:] = flow.y
        self.t = self.data.createVariable('t','float32',('t'))

        # set up the vorticity
        self.w = self.data.createVariable('w','float32',('t','y','x'))
        self.w.setncattr('units','1/s')

        # counter
        self.c = 0

        # add initial flow data
        self.add(flow)


    def add(self, flow) -> None:

        # set the time
        self.t[self.c] = flow.time

        # write vorticity
        self.w[self.c,:,:] = np.fft.irfft2(flow.wh, axes=(-2,-1))

        # update counter
        self.c += 1


    def close(self) -> None:
        # close the Dataset, not mandatory
        self.data.close()