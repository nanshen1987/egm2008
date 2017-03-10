#include <iostream>
#include <fstream>
#include <cassert>
#include <tuple>
#include <map>
#include <algorithm>
#include <cstring>
#include <string>
#include <stdexcept>

/* 
 * const variables; these are pre-set for all three grid files, i.e. 
 * Dg01_cnt2.5x2.5_EGM08_to2190_WGS84_ell_nh
 * eta_cnt2.5x2.5_EGM08_to2190_WGS84_ell_nh
 * xi_cnt2.5x2.5_EGM08_to2190_WGS84_ell_nh
 *
 * The values are extracted from the README file:
 * http://earth-info.nga.mil/GandG/wgs84/gravitymod/egm2008/Anomaly_DOV/README_FIRST_new.pdf
 * and are vital (used throughout) for the program.
 */

constexpr
double beg_lat{90.e0}, end_lat{-90.e0}, //  default latitide limits,
                                        //+ from 90 to -90 (degrees)
       beg_lon{0.e0},  end_lon{360.e0}; //  default longtitude limits,
constexpr
double dlatg  {2.5e0/60.e0},            //  latitude step in degrees
       dlong  {2.5e0/60.e0};            //  longtitude step in degrees
constexpr
int    nrowsg {4320},                   //  # of rows (i.e. latitude point
                                        //+ values)
       ncolsg {8640};                   //  # of cols (i.e. longtitude point
                                        //+ values)
constexpr
unsigned int row_bt_size {sizeof(float)*ncolsg+2*sizeof(int)};
                                        //  size of each row in bytes
    
//  check the limits of the specified (or default) grid
static_assert( beg_lat > end_lat, "" );
static_assert( beg_lon < end_lon, "" );

typedef std::tuple<double, double>         dpair;
typedef std::tuple<long, long>             lpair;
typedef std::map<std::string, std::string> str_str_map;

//  Given a latitude range, return the (inclusive) row numbers they
//+ correspond to.
lpair latitude_line(const dpair&);
//  Given a longtitude range, return the (inclusive) row numbers they
//+ correspond to.
lpair longtitude_line(const dpair&);
//  Parse command line arguments to a string-string map
int cmd_parse(int, char* argc[], str_str_map&);
//  Parse a region of type 'lat0/lat1/lon0/lon1' (i.e. a string)
int limit_parser(const std::string&, double&, double&, double&, double&);
//  perform bilinear interpolation
double
bilinear_interpolation(const double lat, const double lon, 
    const double lat1, const double lat2, 
    const double lon1, const double lon2,
    const double q12, const double q22, const double q11, const double q21);
// help message, usage and epilog
void help();
void usage();
void epilog();

//  for easy-use
struct Point { double lat, lon, val; };

int main(int argc, char* argv[])
{
    //  a dictionary with (any) default options
    str_str_map arg_dict;
    arg_dict["region"]              = std::string( "90.0e0/-90.0e0/0.0e0/360.0e0" );
    arg_dict["point_interpolation"] = std::string("n");

    //  Parse command line arguments to the dictionary
    int status = cmd_parse(argc, argv, arg_dict);
    if (status > 0) {
        std::cerr<<"\n[ERROR] Invalid command line arguments.\n";
        return 1;
    } else if (status < 0) {
        return 0;
    }

    //  open the input file; exit with error if something goes wrong
    auto it = arg_dict.find("input_file");
    if ( it == arg_dict.end() )
    {
        std::cerr<<"\n[ERROR] Must provide an input file!\n";
        return 1;
    }
    std::ifstream fin (it->second.c_str(), std::ios::binary);
    if ( !fin.is_open() )
    {
        std::string fn (argv[1]);
        std::cerr<<"\n[ERROR] Could not open input file \""<< fn <<"\"\n.";
        return 1;
    }

    //  see where we are reporting results; if neccesary open an output file.
    //  From here on, report to 'fp'.
    it = arg_dict.find("output_file");
    std::ostream* fp = &std::cout;
    std::ofstream fout;
    if ( it != arg_dict.end() )
    {
        fout.open(it->second.c_str());
        fp = &fout;
    }

    float  /*exclud {9999.e0},*/ val;
    double rlat, //  latitude of current row
           clon; //  longtitude of current row
#ifdef DEBUG
    double minval {9999.0e0}, minlat, minlon,
           maxval {-9999.0e0},maxlat, maxlon;
    unsigned long values_read {0};
#endif

#ifdef DEBUG
    //  check the file size and report
    fin.seekg(0, std::ios::end);
    auto file_size = fin.tellg();
    std::cout<<"\n** Size of file is "<<file_size;
    std::cout<<"\n** This should be equal to:"
    <<"\n\trows * cols * sizeof(float) + rows * 2 * sizeof(int) ="
    <<nrowsg<<" * "<<ncolsg<<" * "<<sizeof(float)<<" + "
    <<nrowsg*ncolsg*sizeof(float)
    <<"\n\t"<<nrowsg<<" * 2 * "<<sizeof(int)<<" = "
    <<nrowsg*ncolsg*sizeof(float)+nrowsg*2*sizeof(int);
#endif

    //  the size of the file should be equal to:
    //+ rows * cols * sizeof(float) + rows * 2 * sizeof(int)
    //+ that is because (due to the files being written in FORTRAN), in each row
    //+ there is an integer denoting the number of data records, both at the
    //+ begining and at the end of the row.
    fin.seekg(0, std::ios::end);
    assert( fin.tellg() == nrowsg*ncolsg*sizeof(float)+nrowsg*2*sizeof(int) );

    //  get the region
    double st_lat{0.0}, e_lat{0.0}, st_lon{0.0}, e_lon{0.0};
    it = arg_dict.find("region");
    if ( limit_parser(it->second, st_lat, e_lat, st_lon, e_lon) )
    {
        std::cerr<<"\n[ERROR] Cannot parse region limits from string \"" << it->second << "\".\n";
        return 1;
    }
    lpair  lat_lines, // starting and ending rows for given latitude region
           lon_lines; // starting and ending cols for given longtitude region

    //  lines and columns that we are interested in (should be read); these
    //+ functions may throw!
    try {
//#if __GNUC__ > 5
        lat_lines = latitude_line( {st_lat, e_lat} );
        lon_lines = longtitude_line( {st_lon, e_lon} );
//#else

//#endif
    } catch (std::runtime_error& e) {
        return 1;
    }

    //  memory buffer for unwanted rows
    char data_row[row_bt_size];
    //  allocate memory buffer to read in the unwanted columns
    char* data_col;
    std::size_t 
                //  size of leading unwanted lines
                st_col_sz {std::get<0>(lon_lines)*sizeof(float)+sizeof(int)},
                //  ending unwanted cols
                e_col_sz  {(ncolsg-std::get<1>(lon_lines)-1)*sizeof(float)+sizeof(int)};
    //  only allocate one block, with max capacity
    std::size_t col_sz = (st_col_sz>e_col_sz)?st_col_sz:e_col_sz;
    data_col = new char[col_sz];

    //  set cout output format (floating point)
    std::cout.setf(std::ios::fixed, std::ios::floatfield);
    std::cout.precision(10);

    //  we only read to interpolate, not extract values
    //  ----------------------------------------------------------------------
    if (arg_dict["point_interpolation"] == "y") {
        Point* nodes = new Point[4];
        int k {0};
        fin.seekg(0, std::ios::beg);
        for (int i=0; i<std::get<0>(lat_lines); i++) //  read and ignore first rows
        {                                            //+ outside the wanted region
            fin.read(data_row, row_bt_size);         
        }
        for (int i=std::get<0>(lat_lines); i<=std::get<1>(lat_lines); i++)
        {
            fin.read(data_col, st_col_sz); //  read in the unwanted colunms
            rlat = 90.e0 - i*dlatg - dlatg*0.5e0;  // compute latitude of row (degrees)
            for (int j=std::get<0>(lon_lines); j<=std::get<1>(lon_lines); j++)
            {
                clon = j*dlong + dlong*0.5e0; // compute longtitude (degrees)
                fin.read((char*)&val, sizeof(float));
#ifdef DEBUG
                *fp << "\n(node): " << rlat << " " << clon << " " << val;
#endif
                nodes[k].lat = rlat;
                nodes[k].lon = clon;
                nodes[k].val = val;
                ++k;
            }
            fin.read(data_col, e_col_sz); //  read ending unwanted cols
        }
        assert(k==4);
        double lon1{nodes[0].lon}, lon2{nodes[1].lon}, lat1{nodes[2].lat},
               lat2{nodes[0].lat}, q12{nodes[0].val}, q22{nodes[1].val},
               q11{nodes[2].val}, q21{nodes[3].val};
        double x = bilinear_interpolation(st_lat, st_lon, lat1, lat2, lon1, 
            lon2, q12, q22, q11, q21);
        delete[] nodes;
        *fp << "\n" << st_lat << " " << st_lon << " " << x;
    }
    // ------------------------------------------------------------------------
    // interpolation done!

    //  we need to extract values, not interpolate
    //  ----------------------------------------------------------------------
    if (arg_dict["point_interpolation"] == "n") {
        //  go to the top of the file and start reading values. Note that at the
        //+ begining and end of every row, the file has an integer record, denoting
        //+ the length of data (this is a fortran thing...)
        fin.seekg(0, std::ios::beg);
        for (int i=0; i<std::get<0>(lat_lines); i++) //  read and ignore first rows
        {                                            //+ outside the wanted region
            fin.read(data_row, row_bt_size);         
        }
        for (int i=std::get<0>(lat_lines); i<=std::get<1>(lat_lines); i++)
        {
            fin.read(data_col, st_col_sz); //  read in the unwanted colunms
            rlat = 90.e0 - i*dlatg - dlatg*0.5e0;  // compute latitude of row (degrees)
            for (int j=std::get<0>(lon_lines); j<=std::get<1>(lon_lines); j++)
            {
                clon = j*dlong + dlong*0.5e0; // compute longtitude (degrees)
                fin.read((char*)&val, sizeof(float));
                *fp << "\n" << rlat << " " << clon << " " << val;
#ifdef DEBUG
                // in debug mode
                ++values_read;
                if (val < minval)
                {
                    minval = val;
                    minlat = rlat;
                    minlon = clon;
                }
                if (val > maxval)
                {
                    maxval = val;
                    maxlat = rlat;
                    maxlon = clon;
                }
#endif
            }
            fin.read(data_col, e_col_sz); //  read ending unwanted cols
        }
    }

    //  dealocate memory
    delete[] data_col;

#ifdef DEBUG
    if (arg_dict["point_interpolation"] == "n") {
        std::cout << "\nInput file:        "<<arg_dict["input_file"];
        std::cout << "\n# of values read:  "<<values_read;
        std::cout << "\nMinimum value:     "<<minval;
        std::cout << "\nLat of min value:  "<<minlat;
        std::cout << "\nLon of min value:  "<<minlon;
        std::cout << "\nMaximim value:     "<<maxval;
        std::cout << "\nLat of max value:  "<<maxlat;
        std::cout << "\nLon of max value:  "<<maxlon;
        std::cout << "\nInput buf position "<<fin.tellg();
    }
#endif

    std::cout<<"\n";
    return 0;
}

//  Given a latitude range, return the (inclusive) row numbers they
//+ correspond to.
//
//  Given latitude limits, return the row numbers which should be read. The
//+ returned region is inclusive (i.e. both the first and the last line should
//+ be read). If the function returns [x1, x2], then we must read all lines
//+ from x1 up to and including x2.
//  Note that the latitude grid (in the original files) are in descending order,
//+ so the first element of the returned tuple will correspond to the largest
//+ given latitude and the second element to the smallest.
//  The input is a tuple of max and min latitude values, which can be given in
//+ any order (they will be re-arranged inside the function if neccessery (thus
//+ the tuple [19.5, 30.5] will return the same result as [30.5, 19.5].
//  The output is a tuple of long integers, corresponding to the input tuple,
//+ so that:
//  given [x1, x2], where x1 > x2, the function will return [y1, y2], where y1
//+ corresponds to the first latitude smaller than x1 and y2 corresponds to the
//+ first latitude larger than x2.
lpair
latitude_line(const dpair& lat_limits)
{
    double slat = std::get<0>(lat_limits);
    double elat = std::get<1>(lat_limits);
    long   islat, ielat;
    double flat, limit;
#ifdef DEBUG
    double slimit, elimit;
#endif

    //  Note that the grid here is in descending order (i.e. [90, -90]), so
    //+ if needed we should re-arrange the order of the given latitude limits
    if (slat < elat) std::swap(slat, elat);
    
    if (slat >= beg_lat) {
        flat  = 0;
        islat = 0;
    } else {
        flat = (90.e0 - dlatg*0.5e0 - slat) / dlatg + 1;
        islat = static_cast<long>(flat);
    }
    limit = 90.e0 - islat*dlatg - dlatg*0.5e0;
    assert( flat >= 0 && slat >= limit && slat-limit<=dlatg );
#ifdef DEBUG
    slimit = limit;
#endif
    
    if (elat > end_lat) {
        flat = (90.e0 - dlatg*0.5e0 - elat) / dlatg;
        ielat = static_cast<long>(flat);
    } else {
        flat  = 0;
        ielat = nrowsg-1;
    }
    limit = 90.e0 - ielat*dlatg - dlatg*0.5e0;
    assert( flat < nrowsg && limit >= elat && limit-elat<=dlatg );
#ifdef DEBUG
    elimit = limit;
#endif

    //  in case the starting and ending limit is the same value (i.e. the input
    //+ pair is one point), then the computed lines will be in opposite order
    //+ (largest line number will be second). So, rearrange!
    if (slat != elat && islat > ielat) {
        std::cerr<<"\n[ERROR] Could not compute latitude lines properly!";
        throw std::runtime_error("latitude_line() error");
    }
    if (islat > ielat) std::swap(islat, ielat);

#ifdef DEBUG
    std::cout<<"\n** Actual latitude region: ["<<elimit<<", "<<slimit
        <<"] given limits: "<<elat<<", "<<slat;
    std::cout<<"\n** Reading lats from line: "<<islat<<" to line: "<<ielat;
#endif

#if __GNUC__ > 5
    return {islat, ielat};
#else
    return std::make_tuple(islat, ielat);
#endif
}

//  Given a longtitude range, return the (inclusive) row numbers they
//+ correspond to.
//
//  Given longtitude limits, return the col numbers which should be read. The
//+ returned region is inclusive (i.e. both the first and the last column should
//+ be read). If the function returns [x1, x2], then we must read all columns
//+ from x1 up to and including x2.
//  The input is a tuple of max and min longtitude values, which can be given in
//+ any order (they will be re-arranged inside the function if neccessery (thus
//+ the tuple [19.5, 30.5] will return the same result as [30.5, 19.5].
//  The output is a tuple of long integers, corresponding to the input tuple,
//+ so that:
//  given [x1, x2], where x1 < x2, the function will return [y1, y2], where y1
//+ corresponds to the first longtitude larger than x1 and y2 corresponds to the
//+ first longitude smaller than x2.
lpair
longtitude_line(const dpair& lon_limits)
{
    double slon = std::get<0>(lon_limits);
    double elon = std::get<1>(lon_limits);
    long   islon, ielon;
    double flon, limit;
#ifdef DEBUG
    double slimit, elimit;
#endif
    
    if (slon > elon) std::swap(slon, elon);

    if (slon <= beg_lon) {
        flon  = 0;
        islon = 0;
    } else {
        flon = (slon - dlong*0.5e0) / dlong + 1;
        islon = static_cast<long>(flon);
    }
    limit = islon*dlong + dlong*0.5e0;
    assert( flon >= 0 && limit >= slon && limit-slon <= dlong );
#ifdef DEBUG
    slimit = limit;
#endif
    
    if (elon < end_lon) {
        flon = (elon - dlong*0.5e0) / dlong;
        ielon = static_cast<long>(flon);
    } else {
        flon  = 0;
        ielon = ncolsg-1;
    }
    limit = ielon*dlong + dlong*0.5e0;
    assert( flon < ncolsg && elon >= limit && limit-elon<=dlong );
#ifdef DEBUG
    elimit = limit;
#endif

    //  in case the starting and ending limit is the same value (i.e. the input
    //+ pair is one point), then the computed lines will be in opposite order
    //+ (largest line number will be second). So, rearrange!
    if (slon != elon && islon > ielon) {
        std::cerr<<"\n[ERROR] Could not compute longtitude lines properly!";
        throw std::runtime_error("longtitude_line() error");
    }
    if (islon > ielon) std::swap(islon, ielon);

#ifdef DEBUG
    std::cout<<"\n** Actual longtitude region: ["<<slimit<<", "<<elimit
        <<"] given limits: "<<slon<<", "<<elon;
    std::cout<<"\n** Reading lons from line: "<<islon<<" to line: "<<ielon;
#endif

#if __GNUC__ > 5
    return {islon, ielon};
#else
    return std::make_tuple(islon, ielon);
#endif
}

//  Parse command line arguments and assigne them to the given dictionary.
int
cmd_parse(int argv, char* argc[], str_str_map& smap)
{
    if (argv == 1)
    {
            help();
            std::cout << "\n";
            usage();
            std::cout << "\n";
            epilog();
            std::cout << "\n";
            return 0;
    }

    bool r_switch{false}, p_switch{false};

    for (int i = 1; i < argv; i++)
    {
        if ( !std::strcmp(argc[i], "-h") || !std::strcmp(argc[i], "--help"))
        {
            help();
            std::cout << "\n";
            usage();
            std::cout << "\n";
            epilog();
            std::cout << "\n";
            return -1;
        }
        else if ( !std::strcmp(argc[i], "-o") )
        {
            if (i+1 >= argv) { return 1; }
            smap["output_file"] = std::string(argc[i+1]);
            ++i;
        }
        else if ( !std::strcmp(argc[i], "-r") )
        {
            if (i+1 >= argv) { return 1; }
            smap["region"] = std::string(argc[i+1]);
            r_switch = true;
            if (p_switch) {
                std::cerr<<"\n[ERROR] Cannot have both \"-r\" and \"-p\" switches.";
                return 1;
            }
            ++i;
        }
        else if ( !std::strcmp(argc[i], "-p") )
        {
            if (i+1 >= argv) { return 1; }
            std::string token = std::string(argc[i+1]);
            auto j = token.find(",");
            if (j == std::string::npos) { return 1; }
            if (j != token.rfind(",") ) { return 1; }
            auto lat = token.substr(0, j);
            auto lon = token.substr(j+1);
            smap["region"] = (lat+"/"+lat+"/"+lon+"/"+lon);
            smap["point_interpolation"] = std::string{"y"};
            p_switch = true;
            if (r_switch) {
                std::cerr<<"\n[ERROR] Cannot have both \"-r\" and \"-p\" switches.";
                return 1;
            }
            ++i;
        }
        else if ( !std::strcmp(argc[i], "-i") )
        {
            if (i+1 >= argv) { return 1; }
            smap["input_file"] = std::string(argc[i+1]);
            ++i;
        }
        else
        {
            std::cerr<<"\nIrrrelevant cmd: \"" << argc[i] << "\". Skipping!";
        }
    }
    return 0;
}

//  Parse a string of type: 'x1/x2/y1/y2' to a region, where:
//+ x1 -> starting latitude
//+ x2 -> ending latitude
//+ y1 -> starting longtitude
//+ y2 -> ending longtitude
//
//  If anything other than false is returned, an error has occured.
int
limit_parser(const std::string& arg_str, double& lat0, double& lat1,
    double& lon0, double& lon1)
{
    // tokenize the string using '/' as delimeter
    double lary[4] = {-9999, -9999, -9999, -9999};
    int j = 0;
    std::size_t i = 0;
    std::string::size_type pos = arg_str.find_first_of('/', i);
    try
    {
        while (pos != std::string::npos)
        {
            lary[j] = std::stod( arg_str.substr(i, pos-i) );
            ++j;
            i = pos + 1;
            pos = arg_str.find_first_of('/', i);
        }
        lary[j] = std::stod( arg_str.substr(i, pos-i) );
    }
    catch (std::invalid_argument& e)
    {
        return 1;
    }

    // see that all values and no more are assigned.
    if ( j!=3 ) { return 1; }

    lat0 = lary[0];
    lat1 = lary[1];
    lon0 = lary[2];
    lon1 = lary[3];

    return 0;
}

//  q12           q22     ^ latitude axis (y)
//+   o----------o        |
//+   |          |        |
//+   |          |        |
//+   |          |        |
//+   o----------o        +---------> longtitude axis (x)
//+ q11           q21
//+
double
bilinear_interpolation(const double lat, const double lon, 
    const double lat1, const double lat2, 
    const double lon1, const double lon2,
    const double q12, const double q22, const double q11, const double q21)
{
    // x-axis is longtitude
    // y-axis is latitude
    double f_xy1 = ((lon2-lon)/(lon2-lon1))*q11 + ((lon-lon1)/(lon2-lon1))*q21;
    double f_xy2 = ((lon2-lon)/(lon2-lon1))*q12 + ((lon-lon1)/(lon2-lon1))*q22;
    return ((lat2-lat)/(lat2-lat1))*f_xy1 + ((lat-lat1)/(lat2-lat1))*f_xy2;
}

void
help()
{
    std::cout<<"\n";
    std::cout<<"\nProgram egmutil Version 1.0-0\n";
    std::cout<<"\nRead in EGM 2008 related grid files and extract a specific region\n"
    "(given by latitude/longtitude limits in degrees). For more information, see\n"
    "http://earth-info.nga.mil/GandG/wgs84/gravitymod/egm2008/anomalies_dov.html\n"
    "The program can also be used to interpolate at a given point.\n"
    "Note that the input files (as published) are written using BIG ENDIAN binary\n"
    "representation. A conversion to LITTLE INDIAN may be needed depending on the\n"
    "host machine. To perform the conversion, one can use the \"make_SE.f\" utility\n"
    "see http://earth-info.nga.mil/GandG/wgs84/gravitymod/egm2008/endian_convert.html";
    return;
}
void
usage()
{
    std::cout << "\n"
    "Usage:\n"
    " -h or --help\n"
    "\tDisplay (this) help message and exit.\n"
    " -i [input file]\n"
    "\tSpecify the input grid file; this can be any of the Egm2008 (global)\n"
    "\tgrid files (e.g. \"2.5 x 2.5-Minute Free-Air Gravity Anomaly Grid\")\n"
    "\tTo access those files, see http://earth-info.nga.mil/GandG/wgs84/gravitymod/egm2008/anomalies_dov.html\n"
    " -o [output file]\n"
    "\tThis switch is optional. If no ouput file is specified, then the program\n"
    "\twill report to STDOUT. Else, the output is directed to the specified file.\n"
    "\tThe format of the output is \"latitude longtitude value\", where latitude\n"
    "\tand longtitude are in decimal degrees.\n"
    " -r [region]\n"
    "\tSpecify the region to extract from the global grid. The region should be\n"
    "\tgiven as \"min_lat/max_lat/min_lon/max_lon\", where all quantities are\n"
    "\tin decimal degrees.\n"
    " -p [point]\n"
    "\tCompute the value (xi, eta or Dg) at a given point. The point is given\n"
    "\tas: \"latitude,longtitude\" in decimal degrees. The program will use the\n"
    "\tbilinear interpolation algorithm to compute the value at the given point\n"
    "\tbased on the surrounding grid.";

    return;
}

void
epilog()
{
    std::cout << "\n"
    "Copyright 2017 National Technical University of Athens.\n\n"
    "This work is free. You can redistribute it and/or modify it under the\n"
    "terms of the Do What The Fuck You Want To Public License, Version 2,\n"
    "as published by Sam Hocevar. See http://www.wtfpl.net/ for more details.\n"
    "\nSend bugs to: \nxanthos[AT]mail.ntua.gr, \ndemanast[AT]mail.ntua.gr \njorgalan@survey.ntua.gr";
    return;
}
