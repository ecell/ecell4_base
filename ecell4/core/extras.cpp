#include "extras.hpp"
#include <iostream>
#include <sstream>

#if defined(HAVE_BOOST_REGEX)
#include <boost/regex.hpp>
#elif defined(WIN32_MSC)
#include <regex>
#else
#include <regex.h>
#endif /* HAVE_BOOST_REGEX */


namespace ecell4
{

namespace extras
{

#ifdef WITH_HDF5
void save_version_information(H5::CommonFG* root, const std::string& version)
{
    if (version.size() > 32)
    {
        throw IllegalArgument("Version info must be shorter than 32 characters.");
    }
    using namespace H5;
    boost::scoped_ptr<DataSet> dataset(
        new DataSet(root->createDataSet("version", H5::StrType(H5::PredType::C_S1, 32), H5::DataSpace(H5S_SCALAR))));
    dataset->write(version.c_str(), dataset->getDataType());
}

std::string load_version_information(const H5::CommonFG& root)
{
    using namespace H5;
    char buf[32];
    const DataSet dataset(DataSet(root.openDataSet("version")));
    dataset.read(buf, dataset.getDataType());
    return std::string(buf);
}
#endif

std::string load_version_information(const std::string& filename)
{
#ifdef WITH_HDF5
    using namespace H5;
    boost::scoped_ptr<H5::H5File>
        fin(new H5::H5File(filename.c_str(), H5F_ACC_RDONLY));
    return load_version_information(*fin);
#else
    return "";
#endif
}

int mystoi(const std::string& s)
{
    std::stringstream ss;
    ss << s;
    int retval;
    ss >> retval;
    return retval;
}

VersionInformation parse_version_information(const std::string& version)
{
#if defined(HAVE_BOOST_REGEX) || defined(WIN32_MSC)
#if defined(HAVE_BOOST_REGEX)
    using namespace boost;
#else /* WIN32_MSC */
    using namespace std::tr1;
#endif /* HAVE_BOOST_REGEX */
    regex reg("^([^-\\.]+-[^-\\.]+-)([0123456789]+)\\.([0123456789]+)\\.([0123456789]+)$");
    smatch result;
    if (!regex_match(version, result, reg))
    {
        throw std::invalid_argument(
            "a wrong version information was given [" + version + "]"); //XXX:
    }

    const std::string header = result.str(1);
    const int majorno = mystoi(result.str(2));
    const int minorno = mystoi(result.str(3));
    const int patchno = mystoi(result.str(4));

    return VersionInformation(header, majorno, minorno, patchno);
#else /* regex.h */
    regex_t reg;
    int errcode = regcomp(
        &reg, "^([^-\\.]+-[^-\\.]+-)([0123456789]+)\\.([0123456789]+)\\.([0123456789]+)$",
        REG_EXTENDED);
    if (errcode != 0)
    {
        char errbuf[100];
        regerror(errcode, &reg, errbuf, sizeof(errbuf));
        regfree(&reg);
        std::cout << "regcompile error: " << errbuf << std::endl;
        throw IllegalState("regcompile error.");
    }

    regmatch_t match[5];
    errcode = regexec(&reg, version.c_str(), 5, match, 0);
    if (errcode != 0)
    {
        char errbuf[100];
        regerror(errcode, &reg, errbuf, sizeof(errbuf));
        regfree(&reg);
        std::cout << "regexec error: " << errbuf << std::endl;
        throw IllegalState("regexec error.");
    }

    const std::string header = version.substr(match[1].rm_so, match[1].rm_eo - match[1].rm_so);
    const int majorno = mystoi(version.substr(match[2].rm_so, match[2].rm_eo - match[2].rm_so));
    const int minorno = mystoi(version.substr(match[3].rm_so, match[3].rm_eo - match[3].rm_so));
    const int patchno = mystoi(version.substr(match[4].rm_so, match[4].rm_eo - match[4].rm_so));

    regfree(&reg);
    return VersionInformation(header, majorno, minorno, patchno);
#endif /* HAVE_BOOST_REGEX */
}

bool check_version_information(const std::string& version, const std::string& required)
{
    const VersionInformation vinfo1(parse_version_information(version));
    const VersionInformation vinfo2(parse_version_information(required));
    return (vinfo1.header == vinfo2.header
        && vinfo1.majorno >= vinfo2.majorno
        && vinfo1.minorno >= vinfo2.minorno
        && vinfo1.patchno >= vinfo2.patchno);
}

} // extras

} // ecell4
