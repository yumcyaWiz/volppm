"""
    This file is part of the Okapi project

    Copyright 2022 Julian Amann

    Licensed under the Apache License, Version 2.0 (the "License");
    you may not use this file except in compliance with the License.
    You may obtain a copy of the License at

        http://www.apache.org/licenses/LICENSE-2.0

    Unless required by applicable law or agreed to in writing, software
    distributed under the License is distributed on an "AS IS" BASIS,
    WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
    See the License for the specific language governing permissions and
    limitations under the License.

    SPDX-License-Identifier: Apache-2.0
"""

cc_library(
    name = "tinyobjloader",
    hdrs = ["tiny_obj_loader.h"],
    copts = select({
        "@platforms//os:windows": [],
        "//conditions:default": ["-Wno-maybe-uninitialized"],
    }),
    visibility = ["//visibility:public"],
)