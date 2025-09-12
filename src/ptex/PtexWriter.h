#ifndef PtexWriter_h
#define PtexWriter_h

/*
PTEX SOFTWARE
Copyright 2014 Disney Enterprises, Inc.  All rights reserved

Redistribution and use in source and binary forms, with or without
modification, are permitted provided that the following conditions are
met:

  * Redistributions of source code must retain the above copyright
    notice, this list of conditions and the following disclaimer.

  * Redistributions in binary form must reproduce the above copyright
    notice, this list of conditions and the following disclaimer in
    the documentation and/or other materials provided with the
    distribution.

  * The names "Disney", "Walt Disney Pictures", "Walt Disney Animation
    Studios" or the names of its contributors may NOT be used to
    endorse or promote products derived from this software without
    specific prior written permission from Walt Disney Pictures.

Disclaimer: THIS SOFTWARE IS PROVIDED BY WALT DISNEY PICTURES AND
CONTRIBUTORS "AS IS" AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING,
BUT NOT LIMITED TO, THE IMPLIED WARRANTIES OF MERCHANTABILITY, FITNESS
FOR A PARTICULAR PURPOSE, NONINFRINGEMENT AND TITLE ARE DISCLAIMED.
IN NO EVENT SHALL WALT DISNEY PICTURES, THE COPYRIGHT HOLDER OR
CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL,
EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO,
PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR
PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND BASED ON ANY
THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT
(INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE
OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGES.
*/

#include "PtexPlatform.h"
#include <map>
#include <vector>
#include <stdio.h>
#include "Ptexture.h"
#include "PtexIO.h"
#include "PtexReader.h"

struct libdeflate_compressor;

PTEX_NAMESPACE_BEGIN

class PtexMainWriter : public PtexWriter {
public:
    PtexMainWriter(const char* path, PtexTexture* tex,
                   Ptex::MeshType mt, Ptex::DataType dt,
                   int nchannels, int alphachan, int nfaces, bool genmipmaps);

    virtual void release();
    virtual bool close(Ptex::String& error);
    virtual bool writeFace(int faceid, const FaceInfo& f, const void* data, int stride);
    virtual bool writeConstantFace(int faceid, const FaceInfo& f, const void* data);

    virtual void setBorderModes(Ptex::BorderMode uBorderMode, Ptex::BorderMode vBorderMode)
    {
        _extheader.ubordermode = uBorderMode;
        _extheader.vbordermode = vBorderMode;
    }
    virtual void setEdgeFilterMode(Ptex::EdgeFilterMode edgeFilterMode)
    {
            _extheader.edgefiltermode = edgeFilterMode;
    }
    virtual void writeMeta(const char* key, const char* value);
    virtual void writeMeta(const char* key, const int8_t* value, int count);
    virtual void writeMeta(const char* key, const int16_t* value, int count);
    virtual void writeMeta(const char* key, const int32_t* value, int count);
    virtual void writeMeta(const char* key, const float* value, int count);
    virtual void writeMeta(const char* key, const double* value, int count);
    virtual void writeMeta(PtexMetaData* data);

    bool ok(Ptex::String& error) {
        if (!_ok) getError(error);
        return _ok;
    }
    void getError(Ptex::String& error) {
        error = (_error + "\nPtex file: " + _path).c_str();
    }

private:
    virtual ~PtexMainWriter();

    DataType datatype() const { return DataType(_header.datatype); }

    struct MetaEntry {
        std::string key;
        MetaDataType datatype;
        std::vector<uint8_t> data;
        MetaEntry() : datatype(MetaDataType(0)), data() {}
    };

    int writeBlank(FILE* fp, int size);
    int writeBlock(FILE* fp, const void* data, int size);
    int addToDataBlock(std::vector<std::byte>& dataBlock, const void* data, int size);
    int writeZipBlock(FILE* fp, const void* data, int size);
    int readBlock(FILE* fp, void* data, int size);
    int copyBlock(FILE* dst, FILE* src, FilePos pos, int size);
    Res calcTileRes(Res faceres);
    void addMetaData(const char* key, MetaDataType t, const void* value, int size);
    void writeConstFaceBlock(FILE* fp, const void* data, FaceDataHeader& fdh);
    void writeFaceBlock(FILE* fp, const void* data, int stride, Res res,
                       FaceDataHeader& fdh);
    void writeFaceData(FILE* fp, const void* data, int stride, Res res,
                       FaceDataHeader& fdh);
    void writeReduction(FILE* fp, const void* data, int stride, Res res);
    void addToMetaDataBlock(std::vector<std::byte>& dataBlock, const MetaEntry& val);
    void setError(const std::string& error) { _error = error; _ok = false; }
    bool storeFaceInfo(int faceid, FaceInfo& dest, const FaceInfo& src, int flags=0);
    void finish();
    void generateReductions();
    void flagConstantNeighorhoods();
    void storeConstValue(int faceid, const void* data, int stride, Res res);
    void writeMetaData(FILE* fp);

    bool _ok;                                // true if no error has occurred
    std::string _error;                      // the error text (if any)
    std::string _path;                       // file path
    std::string _tilepath;                   // temp tile file path ("<path>.tiles.tmp")
    FILE* _tilefp;                           // temp tile file handle
    Header _header;                          // the file header
    ExtHeader _extheader;                    // extended header
    int _pixelSize;                          // size of a pixel in bytes
    std::vector<MetaEntry> _metadata;        // meta data waiting to be written
    std::map<std::string,int> _metamap;      // for preventing duplicate keys
    libdeflate_compressor* _compressor;      // libzip compression stream

    PtexUtils::ReduceFn* _reduceFn;

    std::string _newpath;                 // path to ".new" file
    std::string _tmppath;                 // temp file path ("<path>.tmp")
    FILE* _tmpfp;                         // temp file handle
    bool _genmipmaps;                     // true if mipmaps should be generated
    std::vector<FaceInfo> _faceinfo;      // info about each face
    std::vector<uint8_t> _constdata;      // constant data for each face
    std::vector<uint32_t> _rfaceids;      // faceid reordering for reduction levels
    std::vector<uint32_t> _faceids_r;     // faceid indexed by rfaceid

    static const int MinReductionLog2 =2; // log2(minimum reduction size) - can tune
    struct LevelRec {
        // note: level 0 is ordered by faceid
        //       levels 1+ are reduction levels (half res in both u and v) and
        //       are ordered by rfaceid[faceid].   Also, faces with a minimum
        //       dimension (the smaller of u or v) smaller than MinReductionLog2
        //       are omitted from subsequent levels.
        std::vector<FilePos> pos;         // position of data blocks within _tmp file
        std::vector<FaceDataHeader> fdh;  // face data headers
    };
    std::vector<LevelRec> _levels;        // info about each level
    std::vector<FilePos> _rpos;           // reduction file positions

    PtexReader* _reader;                  // reader for accessing existing data in file
};


PTEX_NAMESPACE_END

#endif
