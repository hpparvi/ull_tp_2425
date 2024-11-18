# Understaning BHUt Algorithm

## Details

Shape of `bodies`, `bodies%p`, `bodies%p%x` is `n_particles`

## Specifying positions of octants

We are dealing with AABBs (Axis-Aligned Bounding Boxes).
This means that we only need to specify the minimum and maximum values of the cornes:
$$r_{\rm{min}} = (x_{\rm{min}}, y_{\rm{min}}, z_{\rm{min}})$$
$$r_{\rm{max}} = (yx{\rm{max}}, y_{\rm{max}}, z_{\rm{max}})$$
And the positions of the vertices are
$$ \mathcal{C} = 
\{ x, y, z | 
x \in \{x_{\rm{min}}, x_{\rm{max}}\} \ \& \  
y \in \{y_{\rm{min}}, y_{\rm{max}}\} \ \& \  
z \in \{z_{\rm{min}}, z_{\rm{max}}\} 
\}$$
