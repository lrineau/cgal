#include <CGAL/ImageIO.h>
#include <iostream>

#define SHOW(attribut) "\n  "#attribut": " << image->attribut
#define SHOWENUM(enumitem) #enumitem"=" << enumitem

int main(int argc, char** argv)
{
  if(argc == 0) return 1;

  _image* image = ::_readImage(argv[1]);
  if(image)
  {
    std::cerr 
      << "Image infos:"
      << "\ndimensions"
      << SHOW(xdim)
      << SHOW(ydim)
      << SHOW(zdim)
      << SHOW(vdim)
      << "\nvoxel size"
      << SHOW(vx)
      << SHOW(vy)
      << SHOW(vz)
      << "\nword size (in bytes)"
      << SHOW(wdim)
      << "\n";
  }
  else
    std::cerr << "\"" << argv[1] << "\" is not a supported file.\n";
  
  if(argc != 9) return 1;
  const int xmin = std::atoi(argv[3]);
  const int xmax = std::atoi(argv[4]);
  const int ymin = std::atoi(argv[5]);
  const int ymax = std::atoi(argv[6]);
  const int zmin = std::atoi(argv[7]);
  const int zmax = std::atoi(argv[8]);

  const int new_n1 = xmax + 1 - xmin;
  const int new_n2 = ymax + 1 - ymin;
  const int new_n3 = zmax + 1 - zmin;
  const int xdim = image->xdim;
  const int ydim = image->ydim;
  const int wdim = image->wdim;

  _image image2 = *image;
  image2.tx = image->vx * xmin;
  image2.ty = image->vy * ymin;
  image2.tz = image->vz * zmin;
  image2.xdim = new_n1;
  image2.ydim = new_n2;
  image2.zdim = new_n3;
  image2.data =  ImageIO_alloc(std::size_t(new_n1)*new_n2*new_n3*wdim);
  const char* const ptr = (const char*)(image->data);
  char* out_ptr =               (char*)(image2.data);
  for(int k = zmin; k <= zmax; ++k) {
    for(int j = ymin; j <= ymax; ++j) {
      const char * const tmp_ptr = ptr +
        wdim * (xmin + (xdim * (j + ydim * std::size_t(k))));
      out_ptr = std::copy(tmp_ptr, tmp_ptr+wdim*new_n1, out_ptr);
    }
  }
  if(out_ptr !=
     ((char*)(image2.data) + (std::size_t(wdim) * new_n1 * new_n2 *new_n3)))
    return 1;
  ::_writeImage(&image2, argv[2]);
  ImageIO_free(image2.data);
  ::_freeImage(image);
  return 0;
}
