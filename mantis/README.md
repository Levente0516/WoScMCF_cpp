void Impl::compute_interception_list() ... Line 1289-1301 were modifed:
added:        std::vector<BoundingBox> original_intercepted_faces_bb = intercepted_faces_bb[v];
and changed                     set(packed.min_x, j, (float)intercepted_faces_bb[v][i * SimdWidth + j].lower.x);
to                     set(packed.min_x, j, (float) original_intercepted_faces_bb[i * SimdWidth + j].lower.x);