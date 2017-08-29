# Benchmarks
- Designed for Xilinx SDSoC 2016.1. 
- The loop tile labled with "tile_comp" is pipelined with II = 1 and can be flattened by HLS. 
- The burst tranfer is inferred by pipelined loop. 


# Note
The HLS tool has some bug when scheduling the loop tile of FSME. The software implementation of the tiled loop has same output as those produced by the orginal loop. Some output data of the hardware implementation is different. However, the communication time is not affected by this mismatch.

