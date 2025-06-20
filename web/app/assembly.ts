// const assembly = {
//   name: "GRCh38",
//   aliases: ["hg38"],
//   sequence: {
//     type: "ReferenceSequenceTrack",
//     trackId: "GRCh38-ReferenceSequenceTrack",
//     adapter: {
//       type: "BgzipFastaAdapter",
//       fastaLocation: {
//         uri: "https://jbrowse.org/genomes/GRCh38/fasta/hg38.prefix.fa.gz",
//       },
//       faiLocation: {
//         uri: "https://jbrowse.org/genomes/GRCh38/fasta/hg38.prefix.fa.gz.fai",
//       },
//       gziLocation: {
//         uri: "https://jbrowse.org/genomes/GRCh38/fasta/hg38.prefix.fa.gz.gzi",
//       },
//     },
//   },
//   refNameAliases: {
//     adapter: {
//       type: "RefNameAliasAdapter",
//       location: {
//         uri: "https://s3.amazonaws.com/jbrowse.org/genomes/GRCh38/hg38_aliases.txt",
//       },
//     },
//   },
// };

const assembly = {
  name: "Fodor",
  aliases: ["Fodor"],
  sequence: {
    type: "ReferenceSequenceTrack",
    trackId: "GRCh38-ReferenceSequenceTrack",
    adapter: {
      type: "IndexedFastaAdapter",
      fastaLocation: {
        uri: "http://localhost:3000/consensus.fa",
      },
      faiLocation: {
        uri: "http://localhost:3000/consensus.fa.fai",
      }
    },
  },
};

export default assembly;
