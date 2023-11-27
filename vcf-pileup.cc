
#define countof(x) \
	(sizeof((x)) / sizeof(*(x)))

#include <cassert>
#include <string>
#include <zlib.h>
#include <tuple>

namespace {

class GZFile {
public:
	explicit GZFile(const char *path, const char *mode = "r")
		: _file{ ::gzopen(path, mode) }
	{}

	~GZFile()
	{
		if (is_open())
			close();
	}

	bool is_open() const
	{ return _file != nullptr; }

	void close()
	{
		assert(is_open());

		std::ignore = ::gzclose(_file);
		_file = nullptr;
	}

	bool is_eof() const
	{
		assert(is_open());
		return ::gzeof(_file) != 0;
	}

	int get()
	{
		assert(is_open());
		return gzgetc(_file);
	}

	bool getline(std::string& line)
	{
		line.clear();

		for (int ch; (ch = gzgetc(_file)) != -1;)
		{
			switch (ch)
			{
				case '\r':
					if ((ch = gzgetc(_file)) != -1)
						;
					else
						return is_eof();
					if (ch == '\n')
						;
					else
						return ::gzungetc(ch, _file) != -1;
					// FALLTHROUGH

				case '\n':
					return true;

				default:
					line.push_back(char(ch));
					break;
			}
		}

		return is_eof() && !line.empty();
	}

private:
	GZFile(const GZFile&) = delete;
	GZFile& operator=(const GZFile&) = delete;

	gzFile _file;

}; // GZFile

} // (anonymous)

#include <htslib/sam.h>

namespace {

template<class Type>
class SAMObject {
public:
	SAMObject()
		: _ptr{nullptr}
	{}

	SAMObject(SAMObject&& other) = default;

	SAMObject(Type ptr)
		: _ptr{ptr}
	{}

	~SAMObject()
	{
		destroy(_ptr);
		_ptr = nullptr;
	}

	void swap(SAMObject& other)
	{ std::swap(_ptr, other._ptr); }

	void reset()
	{ SAMObject().swap(*this); }

	void reset(Type ptr)
	{ SAMObject(ptr).swap(*this); }

	operator Type() const
	{ return _ptr; }

protected:
	void destroy(Type);

private:
	SAMObject(const SAMObject&) = delete;
	SAMObject& operator=(const SAMObject&) = delete;

	Type _ptr;

}; // SAMObject

template<>
void
SAMObject<samFile *>::destroy(samFile *ptr)
{
	if (ptr != nullptr)
		std::ignore = sam_close(ptr);
}

template<>
void
SAMObject<sam_hdr_t *>::destroy(sam_hdr_t *hdr)
{ sam_hdr_destroy(hdr); }

template<>
void
SAMObject<hts_idx_t *>::destroy(hts_idx_t *idx)
{ hts_idx_destroy(idx); }

template<>
void
SAMObject<hts_itr_t *>::destroy(hts_itr_t *itr)
{ hts_itr_destroy(itr); }

class SAMIterator {
	friend class SAMFile;

public:
	SAMIterator()
		: SAMIterator{nullptr}
	{}

	SAMIterator(SAMIterator&& other) = default;

	~SAMIterator() = default;

protected:
	SAMIterator(hts_itr_t *itr)
		: _itr{itr}
	{}

private:
	SAMIterator(const SAMIterator&) = delete;
	SAMIterator& operator=(const SAMIterator&) = delete;

	SAMObject<hts_itr_t *> _itr;

}; // SAMIterator

template<class Unsigned>
static
Unsigned
adjust_qual(Unsigned qual, bool illumina13)
{
        // NB. from mpileup flag -6 / --illumina1.3+
        if (illumina13)
                qual = qual > 31 ? qual - 31 : 0;

        return qual;
}

static
char
qual2asc(std::uint8_t qual)
{ return qual + 33 < 126 ? qual + 33 : 126; }

class SAMRecordView {
public:
	SAMRecordView()
		: SAMRecordView{nullptr}
	{}

	SAMRecordView(const SAMRecordView&) = default;

	SAMRecordView& operator=(const SAMRecordView&) = default;

	~SAMRecordView() = default;

	const char *qname() const
	{ return bam_get_qname(_bam); }

	// NOTE: or string_view{ bam_get_qname(_bam), _bam->core.l_qname-_bam.core.l_extranul } for qname()
	 // this is '\0' terminated so either works

	std::uint16_t flag() const
	{ return _bam->core.flag; }

	// TODO: individual flags

	int rname_tid() const
	{ return _bam->core.tid; }

	hts_pos_t pos() const
	{ return _bam->core.pos; }

	std::uint8_t mapq() const
	{ return _bam->core.qual; }

	// TODO: raw cigar()

	std::string cigar_str() const
	{
		std::string cigar;

		for (std::uint32_t i = 0; i < _bam->core.n_cigar; ++i)
		{
			cigar += std::to_string( bam_cigar_oplen( bam_get_cigar(_bam)[i] ) );
			cigar +=                 bam_cigar_opchr( bam_get_cigar(_bam)[i] );
		}

		return cigar;
	}

	int rnext_tid() const
	{ return _bam->core.mtid; }

	hts_pos_t pnext() const
	{ return _bam->core.mpos; }

	hts_pos_t tlen() const
	{ return _bam->core.isize; }

	// TODO: raw seq()

	std::string seq_str() const
	{
		std::string seq( _bam->core.l_qseq, '\0' );

		for (int i = 0; i < _bam->core.l_qseq; ++i)
			seq[i] = seq_nt16_str[bam_seqi(bam_get_seq(_bam), i)];

		return seq;
	}

	// TODO: raw qual()

	std::string qual_str() const
	{
		std::string qual( _bam->core.l_qseq, '\0' );

		for (int i = 0; i < _bam->core.l_qseq; ++i)
			qual[i] = qual2asc( bam_get_qual(_bam)[i] );

		return qual;
	}

	hts_pos_t _cigar_len() const
	{ return _bam->core.n_cigar; }

	unsigned _cigar_oplen(hts_pos_t i) const
	{ return bam_cigar_oplen( bam_get_cigar(_bam)[i] ); }

	unsigned _cigar_opchr(hts_pos_t i) const
	{ return bam_cigar_opchr( bam_get_cigar(_bam)[i] ); }

	unsigned _cigar_op(hts_pos_t i) const
	{ return bam_cigar_op( bam_get_cigar(_bam)[i] ); }

	hts_pos_t _seq_len() const
	{ return _bam->core.l_qseq; }

	char _seq_chr(hts_pos_t i) const
	{ return seq_nt16_str[bam_seqi( bam_get_seq(_bam), i )]; }

	std::uint8_t _qual_val(hts_pos_t i) const
	{ return bam_get_qual(_bam)[i]; }

	char _qual_chr(hts_pos_t i) const
	{ return qual2asc( bam_get_qual(_bam)[i] ); }

protected:
	SAMRecordView(bam1_t *bam)
		: _bam{bam}
	{}

private:
	bam1_t *_bam;

}; // SAMRecordView

class SAMRecord : public SAMRecordView {
	friend class SAMFile;

public:
	SAMRecord()
		: SAMRecordView{ &_bam }
		, _bam{}
	{ bam_set_mempolicy( &_bam, BAM_USER_OWNS_STRUCT ); }

	~SAMRecord()
	{ bam_destroy1(&_bam); }

private:
	SAMRecord(const SAMRecord&) = delete;
	SAMRecord& operator=(const SAMRecord&) = delete;

	bam1_t _bam;

}; // SAMRecord

class SAMFile {
	friend class SAMRecordView;

public:
	explicit SAMFile(const char *path, const char *mode)
		: _sam{ sam_open(path, mode) }
		, _hdr{}
		, _idx{}
	{
		if (_sam && *mode != 'w')
			_hdr.reset( sam_hdr_read(_sam) );
		if (_sam && *mode != 'w')
			_idx.reset( sam_index_load(_sam, path) );
	}

	~SAMFile() = default;

	bool is_open() const
	{ return _sam != nullptr; }

	SAMIterator querys(const char *qry)
	{ return sam_itr_querys(_idx, _hdr, qry); }

	SAMIterator queryi(int tid, hts_pos_t beg, hts_pos_t end)
	{ return sam_itr_queryi(_idx, tid, beg, end); }

	SAMIterator queryi(const char *ref, hts_pos_t beg, hts_pos_t end)
	{ return queryi(sam_hdr_name2tid(_hdr, ref), beg, end); }

	bool read(SAMRecord& rec)
	{
		assert(is_open());
		return sam_read1(_sam, _hdr, &rec._bam);
	}

	bool read(SAMRecord& rec, SAMIterator& itr)
	{
		assert(is_open());
		return sam_itr_next(_sam, itr._itr, &rec._bam) >= 0;
	}

	const char *tid2name(int tid, int ref_tid = -1) const
	{
		if (0 <= tid && tid < (*_hdr).n_targets)
			return tid == ref_tid ? "=" : sam_hdr_tid2name(_hdr, tid);
		else
			return "*";
	}

private:
	SAMFile(const SAMFile&) = delete;
	SAMFile& operator=(const SAMFile&) = delete;

	SAMObject<samFile *> _sam;
	SAMObject<sam_hdr_t *> _hdr;
	SAMObject<hts_idx_t *> _idx;

}; // SAMFile

} // (anonymous)

#include <algorithm>
#include <cstring>
#include <string>

namespace {

template<class Type>
class array_view {
public:
	array_view()
		: array_view(nullptr, nullptr)
	{}

	array_view(const array_view&) = default;

	array_view(Type *first, Type *last)
		: _first{first}
		, _last{last}
	{}

	array_view& operator=(const array_view&) = default;

	~array_view() = default;

	Type *begin() const
	{ return _first; }

	Type *end() const
	{ return _last; }

	int compare(array_view& other) const
	{
		return std::lexicographical_compare( other.begin(), other.end(), begin(), end() ) -
			std::lexicographical_compare( begin(), end(), other.begin(), other.end() );
	}

private:
	Type *_first, *_last;

}; // array_view

class string_view : public array_view<const char> {
public:
	string_view() = default;

	string_view(const string_view&) = default;

	string_view(const char *first, const char *last)
		: array_view<const char>{ first, last }
	{}

	string_view(const std::string& string)
		: string_view{ string.data(), string.size() }
	{}

	string_view(const char *data)
		: string_view{ data, std::strlen(data) }
	{}

	string_view(const char *data, std::size_t size)
		: string_view{ &data[0], &data[size] }
	{}

	string_view& operator=(const string_view&) = default;

	~string_view() = default;

	std::string str() const
	{ return std::string{ begin(), end() }; }

}; // string_view

} // (anonymous)

#include <string>
#include <vector>

namespace {

template<class Unsigned, class String>
static
bool
parse_unsigned(Unsigned& value, const String& string)
{
	bool empty = true;

	value = 0;

	for (char ch : string)
	{
		if (!('0' <= ch && ch <= '9'))
			return false;

		Unsigned digit = ch - '0';
		if (!( value <= ( Unsigned(-1) - digit ) / 10 ))
			return false;

		value = 10 * value + digit;
		empty = false;
	}

	return !empty;
}

template<class String>
static
std::vector<string_view>
tokenize(std::vector<string_view>& fields, const String& line, char sep)
{
	fields.clear();
	fields.emplace_back( &*line.begin(), &*line.begin() );

	for (const char& ch : line)
		if (ch == sep)
			fields.emplace_back( &(&ch)[1], &(&ch)[1] );
		else
			fields.back() = { fields.back().begin(), &(&ch)[1] };

	return fields;
}

class VCFReader {
public:
	VCFReader()
		: _line{}
		, _fields{}
	{}

	explicit VCFReader(GZFile& file)
		: VCFReader{}
	{
		while (file.getline(_line))
			if (_line.size() > 0 && _line[0] == '#')
			{
				if (_line.size() > 1 && _line[1] == '#')
					;
				else
				{
					static const char *const expect[] = { "CHROM", "POS", "ID", "REF", "ALT", "QUAL", "FILTER", "INFO" };

					tokenize(_fields, string_view{ &*_line.begin() + 1, &*_line.end() }, '\t');
					if (!(_fields.size() >= countof(expect) && std::equal(
							&expect[0], &expect[countof(expect)], _fields.begin(), [](string_view lhs, string_view rhs)
							{ return lhs.compare(rhs) == 0; } )))
					{
						_throw_error("invalid VCF header");
						return;
					}

					break;
				}
			}
			else
			{
				_throw_error("missing VCF header");
				return;
			}
	}

	~VCFReader() = default;

	const char *error() const
	{ return _fields.empty() ? _line.c_str() : nullptr; }

	bool read(GZFile& file)
	{
		if (error())
			return false;

		while (file.getline(_line))
		{
			std::size_t width = _fields.size();

			tokenize(_fields, _line, '\t');
			if (_fields.size() < width)
				return _throw_error("missing VCF fields");
			if (_fields.size() > width)
				return _throw_error("excess VCF fields");

			return true;
		}

		if (!file.is_eof())
			return _throw_error("read error");

		return false;
	}

	string_view chrom() const
	{ return _fields[0]; }

	string_view pos_str() const
	{ return _fields[1]; }

	long pos() const
	{
		unsigned long pos;

		if (parse_unsigned(pos, pos_str()))
			return long( pos - 1 );
		else
			return -1;
	}

	string_view ref() const
	{ return _fields[3]; }

	string_view alt() const
	{ return _fields[4]; }

protected:
	bool _throw_error(const char *error)
	{
		_line = error;
		_fields.clear();

		return false;
	}

private:
	VCFReader(const VCFReader&) = delete;
	VCFReader& operator=(const VCFReader&) = delete;

	std::string _line;
	std::vector<string_view> _fields;

}; // VCFReader

static
bool
equal_nt(unsigned char nt1, unsigned char nt2)
{
	unsigned char v1 = seq_nt16_table[nt1];
	unsigned char v2 = seq_nt16_table[nt2];

	return v1 == v2 || v1 == seq_nt16_table[unsigned('N')] || v2 == seq_nt16_table[unsigned('N')];
}

template<class Iterator1, class Iterator2>
static
bool
equal_nt_seq(Iterator1 it1, Iterator1 end1, Iterator2 it2)
{
	for (; it1 != end1; ++it1, ++it2)
		if (!equal_nt( *it1, *it2 ))
			return false;

	return true;
}

static
string_view
basename(string_view path)
{
	const char *p = &*path.begin();
	const char *q = &*path.end();

	while (p != q && q[-1] == '/')
		--q;
	for (const char *z = p; z != q; ++z)
		if (*z == '/')
			p = &z[1];

	return string_view{ p, q };
}

static
std::pair<string_view, string_view>
splitext(string_view path)
{
	const char *p = path.begin();

	for (const char *z = p; z != path.end(); ++z)
		if (*z == '.')
			p = z;

	return std::make_pair( string_view{ path.begin(), p }, string_view{ p, path.end() } );
}

} // (anonymous)

#include <cstdio>
#include <cstdlib>
#include <list>
#include <map>

int
main(int argc, char **argv)
{
	bool header = false;
	bool illumina13 = false;
	std::uint8_t min_mq = 0;
	std::uint8_t min_bq = 13;

	array_view<const char *const> args{ &argv[1], &argv[argc] };

	for (; args.end() - args.begin() >= 1 && ( *args.begin() )[0] == '-'; args = { args.begin() + 1, args.end() })
		     if (string_view( *args.begin() ).str() == "-h")
			header = true;
		else if (string_view( *args.begin() ).str() == "-6")
			illumina13 = true;
		else if (string_view( *args.begin() ).str() == "-q")
		{
			if (parse_unsigned( min_mq, string_view(*( args.begin() + 1 )) ))
				args = { args.begin() + 1, args.end() };
			else
				goto usage;
		}
		else if (string_view( *args.begin() ).str() == "-Q")
		{
			if (parse_unsigned( min_bq, string_view(*( args.begin() + 1 )) ))
				args = { args.begin() + 1, args.end() };
			else
				goto usage;
		}
		else
			goto usage;

	if (args.end() - args.begin() < 2)
	{
usage:
		std::fprintf(stderr, "Usage: %s [-h] [-6] [-q mq] [-Q bq] input.vcf input.bam [...]" "\n", argv[0]);
		return EXIT_FAILURE;
	}

	const char *vcf_fn = *args.begin();

	array_view<const char *const> bam_fns{ args.begin() + 1, args.end() };
	std::list<SAMFile> bams;

	for (const char *bam_fn : bam_fns)
	{
		bams.emplace_back(bam_fn, "r");
		if (!bams.back().is_open())
		{
			std::fprintf(stderr, "%s: %s: %s" "\n", argv[0], bam_fn, "failed to open for read");
			return EXIT_FAILURE;
		}
	}

	GZFile vcf_fp{ vcf_fn, "r" };
	if (!vcf_fp.is_open())
	{
		std::fprintf(stderr, "%s: %s: %s" "\n", argv[0], vcf_fn, "failed to open for read");
		return EXIT_FAILURE;
	}

	{
		if (header)
		{
			std::printf("##fileformat=VCFv4.1\n");

			std::printf("##FORMAT=<ID=AD,Number=R,Type=Integer,Description=\"Read depth uniquely matching each allele\">\n");
			std::printf("##FORMAT=<ID=ADF,Number=R,Type=Integer,Description=\"Read depth in the forward strand\">\n");
//			std::printf("##FORMAT=<ID=ADR,Number=R,Type=Integer,Description=\"Read depth in the backward strand\">\n");
			std::printf("##FORMAT=<ID=XD,Number=1,Type=Integer,Description=\"Extra number of reads ambiguously matching alleles\">\n");
			std::printf("##FORMAT=<ID=XDF,Number=1,Type=Integer,Description=\"Extra reads in the forward strand\">\n");
//			std::printf("##FORMAT=<ID=XDR,Number=1,Type=Integer,Description=\"Extra reads in the backward strand\">\n");
			std::printf("##FORMAT=<ID=DP,Number=1,Type=Integer,Description=\"Read depth (quality filtered, regardless of allele)\">\n");
			std::printf("##FORMAT=<ID=DPF,Number=1,Type=Integer,Description=\"Read depth in the forward strand\">\n");
//			std::printf("##FORMAT=<ID=DPR,Number=1,Type=Integer,Description=\"Read depth in the backward strand\">\n");
			std::printf("##FORMAT=<ID=RDP,Number=1,Type=Integer,Description=\"Raw read depth (unfiltered)\">\n");

			// TODO: print program version/date
		}

		std::printf("#" "%s" "\t" "%s" "\t" "%s" "\t" "%s" "\t" "%s" "\t" "%s" "\t" "%s" "\t" "%s" "\t" "%s",
			"CHROM", "POS", "ID", "REF", "ALT", "QUAL", "FILTER", "INFO", "FORMAT");

		for (const char *bam_fn : bam_fns)
			std::printf("\t" "%s", splitext(basename(bam_fn)).first.str().c_str());

		std::printf("\n");
	}

	VCFReader vcf_reader{ vcf_fp };
	while (vcf_reader.read( vcf_fp ))
	{
		string_view chrom = vcf_reader.chrom();
		hts_pos_t pos = vcf_reader.pos();
		string_view ref = vcf_reader.ref();
		string_view alt = vcf_reader.alt();

		std::printf("%.*s" "\t" "%ld" "\t" "." "\t" "%.*s" "\t" "%.*s" "\t" "." "\t" "." "\t" "." "\t" "%s",
			int( chrom.end() - chrom.begin() ), chrom.begin(),
			long( pos + 1 ),
			int( ref.end() - ref.begin() ), ref.begin(),
			int( alt.end() - alt.begin() ), alt.begin(),
			"AD" ":" "ADF" /* ":" "ADR" */ ":" "XD" ":" "XDF" /* ":" "XDR" */ ":" "DP" ":" "DPF" /* ":" "DPR" */ ":" "RDP");

		for (SAMFile& bam : bams)
		{
			SAMIterator itr = bam.queryi( chrom.str().c_str(), pos, pos + ( ref.end() - ref.begin() ) );

			std::size_t cnt[2][2][2] = { { { 0 } } };
			std::size_t cnt_rej = 0;

			for (SAMRecord rec; bam.read(rec, itr);)
				if (rec.mapq() >= min_mq)
				{
					auto consume_ref = [](unsigned op) { return ( bam_cigar_type(op) & 0x2 ) != 0; };
					auto consume_qry = [](unsigned op) { return ( bam_cigar_type(op) & 0x1 ) != 0; };

					std::string qry;
					bool match_ref;
					bool match_alt;
					{
					hts_pos_t p = rec.pos();
					hts_pos_t q = 0;

					bool bind_l = false;
					bool bind_r = false;

					for (hts_pos_t i = 0; i < rec._cigar_len(); ++i)
						for (hts_pos_t l = 0; l < rec._cigar_oplen(i); ++l)
						{
							bool con_ref = consume_ref( rec._cigar_op(i) );
							bool con_qry = consume_qry( rec._cigar_op(i) );

							bool is_clip = rec._cigar_op(i) == BAM_CREF_SKIP || rec._cigar_op(i) == BAM_CSOFT_CLIP;

							if (p < pos)
							{
								if (!is_clip)
									bind_l = true;
							}
							else if (p < pos + ( ref.end() - ref.begin() ) + !( con_ref == con_qry ))
							{
								if (con_qry && !is_clip)
								{
									if (adjust_qual( rec._qual_val(q), illumina13 ) >= min_bq)
										qry.push_back( rec._seq_chr(q) );
									else
										qry.push_back( 'N' );
								}
							}
							else
							{
								if (!is_clip)
									bind_r = true;
							}

							if (con_ref)
								++p;
							if (con_qry)
								++q;
						}

					auto match = [](string_view lhs, string_view rhs, bool bind_l, bool bind_r) {
						auto size = [](string_view str) {
							return std::size_t( std::end(str) - std::begin(str) );
						};

						if (bind_l && bind_r)
							return size(lhs) == size(rhs) && equal_nt_seq( lhs.begin(), lhs.end(), rhs.begin() );
						else if (bind_l)
							return size(lhs) >= size(rhs) && equal_nt_seq( lhs.begin(), lhs.begin() + size(rhs), rhs.begin() );
						else if (bind_r)
							return size(lhs) >= size(rhs) && equal_nt_seq( lhs.end() - size(rhs), lhs.end(), rhs.begin() );
						else
						{
							// TODO: O(n+m) please?
							for (std::size_t s = 0; s + size(rhs) <= size(lhs); ++s)
								if (equal_nt_seq( lhs.begin() + s, lhs.begin() + s + size(rhs), rhs.begin() ))
									return true;

							return false;
						}
					};

					match_ref = match( ref, qry, bind_l, bind_r );
					match_alt = match( alt, qry, bind_l, bind_r );
					}

					++cnt[match_alt][match_ref][ ( rec.flag() & BAM_FREVERSE ) != 0 ];
				}
				else
					++cnt_rej;

			std::printf("\t" "%zu,%zu" ":" "%zu,%zu" /* ":" "%zu,%zu" */ ":" "%zu" ":" "%zu" /* ":" "%zu" */ ":" "%zu" ":" "%zu" /* ":" "%zu" */ ":" "%zu",
				                              cnt[0][1][0] + cnt[0][1][1]                                                                       , // AD REF
				                                                            cnt[1][0][0] + cnt[1][0][1]                                         , // AD ALT
				                              cnt[0][1][0]                                                                                      , // ADF REF
				                                                            cnt[1][0][0]                                                        , // ADF ALT
//				                                             cnt[0][1][1]                                                                       , // ADR REF
//				                                                                           cnt[1][0][1]                                         , // ADR ALT
				cnt[0][0][0] + cnt[0][0][1]                                                                                                     , // XD
				cnt[0][0][0]                                                                                                                    , // XDF
//				               cnt[0][0][1]                                                                                                     , // XDR
				cnt[0][0][0] + cnt[0][0][1] + cnt[0][1][0] + cnt[0][1][1] + cnt[1][0][0] + cnt[1][0][1] + cnt[1][1][0] + cnt[1][1][1]           , // DP
				cnt[0][0][0] +                cnt[0][1][0] +                cnt[1][0][0] +                cnt[1][1][0]                          , // DPF
//				               cnt[0][0][1] +                cnt[0][1][1] +                cnt[1][0][1] +                cnt[1][1][1]           , // DPR
				cnt[0][0][0] + cnt[0][0][1] + cnt[0][1][0] + cnt[0][1][1] + cnt[1][0][0] + cnt[1][0][1] + cnt[1][1][0] + cnt[1][1][1] + cnt_rej ); // RDP
		}

		std::printf("\n");
	}

	if (!vcf_fp.is_eof())
	{
		std::fprintf(stderr, "%s: %s: %s" "\n", argv[0], vcf_fn, vcf_reader.error());
		return EXIT_FAILURE;
	}

	return EXIT_SUCCESS;
}
