#include "../headers/read.hpp"


// unnamed namespace
// don't really want these functions floating around
namespace {

  // trim helpers
  inline void ltrim(std::string& s) {
    size_t i = 0;
    while (i < s.size() && std::isspace(static_cast<unsigned char>(s[i]))) ++i;
    s.erase(0, i);
  }

  inline void rtrim(std::string& s) {
    size_t i = s.size();
    while (i > 0 && std::isspace(static_cast<unsigned char>(s[i - 1]))) --i;
    s.erase(i);
  }

  inline std::string trim_copy(std::string s) {
    ltrim(s);
    rtrim(s);
    return s;
  }

  // reads next meaningful value line: strips # comments, trims, skips empties
  class LineReader {
  public:
    explicit LineReader(const std::string& path) : in_(path) {
      if (!in_) throw std::runtime_error("Could not open config: " + path);
    }

    std::string next_value() {
      std::string line;
      while (std::getline(in_, line)) {
        line_no_++;

        // strip comment
        if (auto pos = line.find('#'); pos != std::string::npos) {
          line.erase(pos);
        }

        line = trim_copy(line);
        if (!line.empty()) return line;
      }
      throw std::runtime_error("Unexpected end of config (line " + std::to_string(line_no_) + ")");
    }

    int line_no() const { return line_no_; }

  private:
    std::ifstream in_;
    int line_no_ = 0;
  };

  // robust int parse via from_chars
  int parse_int(const std::string& s, const char* what, int line_no) {
    int v{};
    auto b = s.data();
    auto e = s.data() + s.size();
    auto [ptr, ec] = std::from_chars(b, e, v);
    if (ec != std::errc{} || ptr != e) {
      throw std::runtime_error(std::string("Invalid integer for ") + what +
                              " at line " + std::to_string(line_no) + ": '" + s + "'");
    }
    return v;
  }

  // float parsing: std::stof with full-consumption check
  float parse_float(const std::string& s, const char* what, int line_no) {
    try {
      size_t idx = 0;
      float v = std::stof(s, &idx);
      if (idx != s.size()) throw std::runtime_error("trailing");
      return v;
    } catch (...) {
      throw std::runtime_error(std::string("Invalid float for ") + what +
                              " at line " + std::to_string(line_no) + ": '" + s + "'");
    }
  }

} // namespace

Config loadConfig(const std::string& path_to_config, const std::string& path_to_home) {
  
  Config c;
  c.path = path_to_home;

  // All the rest is read in from config.asc file

  LineReader r(c.path+path_to_config);

  // Must match the exact order of config file
  c.aeroclim = r.next_value();
  c.xadim    = parse_int(r.next_value(), "xadim", r.line_no());
  c.yadim    = parse_int(r.next_value(), "yadim", r.line_no());

  c.hclim = r.next_value();
  c.xhdim = parse_int(r.next_value(), "xhdim", r.line_no());
  c.yhdim = parse_int(r.next_value(), "yhdim", r.line_no());

  c.o3clim = r.next_value();
  c.xo3dim = parse_int(r.next_value(), "xo3dim", r.line_no());
  c.yo3dim = parse_int(r.next_value(), "yo3dim", r.line_no());

  c.landuse_image  = r.next_value();

  c.use_modis_brdf_albedo = parse_int(r.next_value(), "use_modis_brdf_albedo", r.line_no());
  c.modis_brdf_dir = r.next_value();

  c.clut_spec      = r.next_value();
  c.clut_spec_h2o  = r.next_value();
  c.clut_spec_o3   = r.next_value();
  c.lambda_cor     = r.next_value();

  c.cloud_list       = r.next_value();
  c.out_global_path= r.next_value();
  c.out_beam_path  = r.next_value();
  c.out_cal_path   = r.next_value();
  c.out_clear_path  = r.next_value();

  c.gr_alb_file = r.next_value();

  c.ggdim     = parse_int(r.next_value(), "ggdim", r.line_no());
  c.ssadim    = parse_int(r.next_value(), "ssadim", r.line_no());
  c.aoddim    = parse_int(r.next_value(), "aoddim", r.line_no());
  c.h2odim    = parse_int(r.next_value(), "h2odim", r.line_no()); // fixed name
  c.o3dim     = parse_int(r.next_value(), "o3dim", r.line_no());
  c.num_bands  = parse_int(r.next_value(), "bandsdim", r.line_no());
  c.latdim    = parse_int(r.next_value(), "latdim", r.line_no()); // fixed name
  c.londim    = parse_int(r.next_value(), "londim", r.line_no()); // fixed name

  c.latbegin  = parse_float(r.next_value(), "begin of lat", r.line_no());
  c.lonbegin  = parse_float(r.next_value(), "begin of lon", r.line_no());
  c.dxy       = parse_float(r.next_value(), "dxy", r.line_no());
  c.deltalon  = parse_float(r.next_value(), "deltalon", r.line_no());
  c.iconflag  = parse_int(r.next_value(), "iconflag", r.line_no());
  c.iconres   = parse_float(r.next_value(), "iconres", r.line_no()); // fixed value source

  // Idiot checks
  if (c.xadim <= 0 || c.yadim <= 0) throw std::runtime_error("Invalid aerosol dimensions");
  if (c.latdim <= 0 || c.londim <= 0) throw std::runtime_error("Invalid grid dimensions");
  if (c.use_modis_brdf_albedo != 0 && c.use_modis_brdf_albedo != 1) {
    throw std::runtime_error("Invalid MODIS BRDF albedo flag; expected 0 or 1");
  }

  return c;
}