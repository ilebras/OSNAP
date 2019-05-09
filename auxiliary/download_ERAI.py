from ecmwfapi import ECMWFDataServer

server = ECMWFDataServer()

# First try to just download eastward wind stress for August 2014...

server.retrieve({
    "class": "ei",
    "dataset": "interim",
    "date": "2014-08-01/to/2016-08-01",
    "expver": "1",
    "grid": "0.75/0.75",
    "levtype": "sfc",
    "param": "iews/inss/146.128/147.128",
    "step": "24",
    "time": "00:00:00",
    "stream": "oper",
    "type": "fc",
    "area" : "75/-50/55/-20",
    "format": "netcdf",
    "target": "../data/aux_data/ERA_1804/tau_hflux_180411.nc",
    })


dat['inss'].mean(dim='time').plot()


dat['iews'].mean(dim='time').plot()
