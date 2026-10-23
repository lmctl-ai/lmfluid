# Libraries and licenses

MFEM 4.8 is pinned at a01719101027383954b69af1777dc828bf795d62 (official v4.8 tag), BSD-3-Clause. The build uses the serial CPU implementation with exceptions. MPI, METIS, hypre, OpenMP, CUDA/HIP, LAPACK, libCEED, GSLIB, zlib and all optional external packages are disabled. No Eigen fallback or OpenFOAM dependency is present.

The full MFEM LICENSE and NOTICE are embedded, along with bundled Gecko (BSD-3), TinyXML2 (zlib), Picojson (BSD-2) and Zstr (MIT) notice texts. Optional AmgX and upstream testing components are disabled. The scripting SDK and adapted expression compiler are MIT and the SDK is pinned by commit and SHA-256 in sdk.lock. Our wrappers and D2Q9 are original MIT code.

Deps are fetched by sync-deps.sh into ignored third_party/. No MFEM source is committed. --verify checks the commit and tracked source cleanliness. No source patches are required unless a platform build demonstrates otherwise. Standard compiler runtimes retain their applicable system/runtime-exception terms.

Primary sources: [MFEM license](https://github.com/mfem/mfem/blob/v4.8/LICENSE), [MFEM serial build](https://mfem.org/building/), [MFEM installation options](https://github.com/mfem/mfem/blob/v4.8/INSTALL).
