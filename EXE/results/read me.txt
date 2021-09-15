The results are float32 binary data:
	.air    is the airscan raw data, dim: [col row]
	.offset is the offset scan,      dim: [col row]
	.scan   is the phantom scan,     dim: [col row view]
	.prep   is the -log(scan/air),   dim: [col row view]

If you want to change filename or path of the results, modify the following line in EXE/cfg/Protocol.cfg:
cfg.resultsName = "./results/simulation_test"

