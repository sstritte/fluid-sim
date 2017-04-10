### Parallel Black and White Photo Colorizer
Jennifer Chou (jtchou), Sienna Stritter (sstritte)

## Proposal
### Summary
We are going to optimize an algorithm that uses deep learning to automatically create plausible color versions of black and white photos. 

### Background
Our project is inspired by the work of Richard Zhang, Phillip Isola, and Alexei A. Efros from UC Berkeley. They configured a convolutional neural network (using a deep learning framework called Caffe) that attempts to automatically colorize a black and white photo. We seek to implement our own version of the algorithm by creating a similar neural network setup and parallelizing it.

We will attempt to distribute the network computation across multiple GPU cores, which will speed up calculations but also introduce difficulties with communication across cores, since the overall network proceeds in layers and terminates according to a global loss value. Parallelism can also be used within a layer to perform operations on each pixel.

### The Challenge



### Resources
Our idea is based on this algorithm:
[Colorful Image Colorization](https://arxiv.org/pdf/1603.08511.pdf)
2016. Richard Zhang, Phillip Isola, and Alexei A. Efros.

We will try to implement the colorizer algorithm using a parallelizable neural network based off of:
[Caffe: Convolutional Architecture for Fast Feature Embedding](https://github.com/BVLC/caffe)
2014. Jia, Yangqing and Shelhamer, Evan and Donahue, Jeff and Karayev, Sergey and Long, Jonathan and Girshick, Ross and Guadarrama, Sergio and Darrell, Trevor.

### Goals and Deliverables

### Platform Choice
We will use CUDA.

### Schedule 


 
