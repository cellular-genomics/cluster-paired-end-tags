#define BOOST_SPIRIT_DEBUG

#include <fstream>
#include <cstdio>

#include <boost/fusion/adapted/struct.hpp>
#include <boost/spirit/include/qi.hpp>
#include <boost/iostreams/device/mapped_file.hpp>
#include <boost/utility/string_ref.hpp>
#include <boost/program_options.hpp>

namespace qi = boost::spirit::qi;
using sref   = boost::string_ref;
namespace po = boost::program_options;
using namespace std;

namespace boost { namespace spirit { namespace traits {
    template <typename It>
    struct assign_to_attribute_from_iterators<sref, It, void> {
        static void call(It f, It l, sref& attr) { attr = { f, size_t(std::distance(f,l)) }; }
    };
} } }

struct PET {
    sref chrom1;
    int start1, end1;
    sref chrom2;
    int start2, end2;
    int count;

    bool operator<(const PET& a) const
    {
        if (chrom1 != a.chrom1) return chrom1 < a.chrom1;
        if (start1 != a.start1) return start1 < a.start1;
        if (end1 != a.end1) return end1 < a.end1;
        if (chrom2 != a.chrom2) return chrom2 < a.chrom2;
        if (start2 != a.start2) return start2 < a.start2;
        if (end2 != a.end2) return end2 < a.end2;
        return false;
    }

    friend ostream& operator<<(ostream& os, const PET& pet)
    {
        os << pet.chrom1 << ":" << pet.start1 << "-" << pet.end1 << " " << pet.chrom2 << ":" << pet.start2 << "-" << pet.end2 << " (" << pet.count << ")";
        return os;
    }    
};

BOOST_FUSION_ADAPT_STRUCT(PET, (sref,chrom1)(int,start1)(int,end1)(sref,chrom2)(int,start2)(int,end2)(int,count))

int main(int ac, char** av) {
    string input_filename, output_filename;
    int pet_cutoff, cluster_cutoff, extension, self_ligation;

    try
    {
        po::options_description desc("Allowed options");
        desc.add_options()
            ("help,h", "help message")
            ("input,i", po::value<string>(&input_filename)->required(), ".bedpe file containing raw PETs (no header)")
            ("output,o", po::value<string>(&output_filename)->required(), ".bedpe output file with clustered PETs (no header)")
            ("self_ligation,s", po::value<int>(&self_ligation)->default_value(8000), "Self-ligation genomic span (default: 8000)")
            ("extension,e", po::value<int>(&extension)->default_value(500), "No of base pairs to extend both ends of PETs (default: 500)")
            ("pet_cutoff,p", po::value<int>(&pet_cutoff)->default_value(2), "Minimum number of PET counts to take PET into consideration (default: 2)")
            ("cluster_cutoff,c", po::value<int>(&cluster_cutoff)->default_value(4), "Minimum number of total counts to consider as a cluster (default: 4)");

        po::variables_map vm;
        po::store(po::parse_command_line(ac, av, desc), vm);
        if (vm.count("help")) {
            cout << desc << "\n";
            return 0;
        }
        po::notify(vm);
    }catch (const po::error &ex)
    {
        cerr << ex.what() << '\n';
        return 0;
    }

    // reading the file
    cout << "Reading PETs from " << input_filename << endl;
    boost::iostreams::mapped_file_source mmap(input_filename);

    // parsing
    using namespace qi;
    vector<PET> pets;
    pets.reserve(44455464);
    if (parse(mmap.begin(), mmap.end(), 
                (raw[*~char_("\t")] >> '\t' >> int_ >> '\t' >> int_ >> '\t' >> raw[*~char_("\t")] >> '\t' >> int_ >> '\t' >> int_ >> '\t' >> int_) % eol, 
                pets))
    {
        cout << "Parsed " << pets.size() << " PETs" << endl;
    } else {
        cout << "Failed after " << pets.size() << " PETs" << endl;
    };

    // sorting
    cout << "Sorting... " << flush;
    sort(pets.begin(), pets.end());
    cout << "Done." << endl;

    // pre-proccess
    cout << "Preprocessing (Extension: "<< extension << "bp, Self-ligation genomic span: " << self_ligation << "bp, PET cutoff: " << pet_cutoff << ")... " << flush;
    for (int i = 0; i<pets.size(); i++){
        PET &pet = pets[i];
        // check the data integrity
        if (pet.chrom1 != pet.chrom2 || pet.start1 > pet.end1 || pet.start2 > pet.end2){
            cout << "Inter-chromosomal or misordered PETs are ignored: #" << (i+1) << " " << pet << endl;
            pet.count = 0;
        }

        // remove self ligating PETs
        if (pet.start2 - pet.end1 < self_ligation)
            pet.count = 0;

        // remove PETs with count < PET cutoff
        if (pet.count < pet_cutoff)
            pet.count = 0;

        // add extension
        pet.start1 = max(0, pet.start1 - extension);
        pet.end1 = pet.end1 + extension;
        pet.start2 = max(0, pet.start2 - extension);
        pet.end2 = pet.end2 + extension;
    }
    cout << "Done." << endl;

    // clustering
    cout << "Clustering... " << flush;
    for (int i =0; i<pets.size(); i++) {
        if (pets[i].count == 0) 
            continue;
        int j = i + 1;
        PET &cluster = pets[i];
        while (j < pets.size()) {
            PET &can = pets[j];
            if (can.start1 > cluster.end1 or can.chrom1 != cluster.chrom1)
                break;
            if
            ((((pets[i].start1 <= pets[j].start1) && (pets[j].start1 <= pets[i].end1)) || 
                 ((pets[i].start1 <= pets[j].end1) && (pets[j].end1 <= pets[i].end1))) &&
                (((pets[i].start2 <= pets[j].start2) && (pets[j].start2 <= pets[i].end2)) || 
                 ((pets[i].start2 <= pets[j].end2) && (pets[j].end2 <= pets[i].end2)))) 
                 {
                cluster.start1 = min(cluster.start1, can.start1);
                cluster.end1 = max(cluster.end1, can.end1);
                cluster.start2 = min(cluster.start2, can.start2);
                cluster.end2 = max(cluster.end2, can.end2);
                cluster.count += can.count;
                can.count = 0;
            }
            j++;            
        }
    }
    cout << "Done." << endl;

    // save to file
    cout << "Saving to " << output_filename << " (cluster cufoff: " << cluster_cutoff << ")... " << flush;
    ofstream ofs(output_filename, ofstream::out);
    char buff[64];
    int c = 0;
    for (int i =0; i<pets.size(); i++){
        PET &cluster = pets[i];
        if (cluster.count >= cluster_cutoff){
            std::string chrom1(cluster.chrom1.data(), cluster.chrom1.size()), chrom2(cluster.chrom2.data(), cluster.chrom2.size());
            sprintf(buff, "%s\t%d\t%d\t%s\t%d\t%d\t%d\n", chrom1.c_str(), cluster.start1, cluster.end1, chrom2.c_str(), cluster.start2, cluster.end2, cluster.count);
            ofs << buff;
            c++;
        }
    }
    ofs.close();
    cout << "Done. Saved " << c << " clusters." << endl;
}