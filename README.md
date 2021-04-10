# Path-tracer
This project is about building a **ray-tracer**, renderering images (exported in ppm format) of 3D scenes based on global illumination, simulating diffuse, specular and dielectric materials, by giving a numerical solution to the rendering equation using **Monte-Carlo** path tracing with **Russian Roulette**.
<br/>Currently, renderer features :<br/>
-specular, diffuse and refraction<br/>
-soft shadows<br/>
-antialiasing (with tent filter)<br/>
-Multi-threaded rendering with OpenMP.<br/><br/>

+to do:<br/>
-implementing anisotropic reflection.<br/>
-handling motion blur.<br/>
-adding depth of field and caustic effects.<br/>
-surface textures.<br/>
-and many more..<br/>
Done:<br/>
-~metallic reflection~<br/>
-~parallelizing the image rendering~<br/>
-~~add primitives (planes, cubes...)~~<br/>
<br/>
<br/>
<br/>Demo: (128 samples) <br/>

![alt text](https://github.com/Yass-g/Path-tracer/blob/master/RayTrace/exple/01_spheres_128samples.png)










