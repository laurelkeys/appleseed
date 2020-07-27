
//
// This source file is part of appleseed.
// Visit https://appleseedhq.net/ for additional information and resources.
//
// This software is released under the MIT license.
//
// Copyright (c) 2020 Tiago Chaves, The appleseedhq Organization
//
// Permission is hereby granted, free of charge, to any person obtaining a copy
// of this software and associated documentation files (the "Software"), to deal
// in the Software without restriction, including without limitation the rights
// to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
// copies of the Software, and to permit persons to whom the Software is
// furnished to do so, subject to the following conditions:
//
// The above copyright notice and this permission notice shall be included in
// all copies or substantial portions of the Software.
//
// THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
// IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
// FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
// AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
// LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
// OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN
// THE SOFTWARE.
//

// Interface header.
#include "tonemappostprocessingstage.h"

// appleseed.renderer headers.
#include "renderer/modeling/frame/frame.h"
#include "renderer/modeling/postprocessingstage/postprocessingstage.h"
#include "renderer/modeling/postprocessingstage/effect/tonemapapplier.h"

// appleseed.foundation headers.
#include "foundation/containers/dictionary.h"
#include "foundation/image/canvasproperties.h"
#include "foundation/image/conversion.h"
#include "foundation/image/image.h"
#include "foundation/math/scalar.h"
#include "foundation/math/vector.h"
#include "foundation/utility/api/specializedapiarrays.h"
#include "foundation/utility/makevector.h"

// Standard headers.
#include <cstddef>
#include <memory>
#include <string>
#include <vector>

using namespace foundation;

namespace renderer
{

namespace
{

    //
    // Tone mapping operators.
    //

    struct ToneMapOperator
    {
        const char* label;
        const char* id;
    };

    constexpr const ToneMapOperator AcesNarkowicz       { "ACES (Narkowicz)",       "aces_narkowicz" };
    constexpr const ToneMapOperator AcesUnreal          { "ACES (Unreal)",          "aces_unreal" };
    constexpr const ToneMapOperator FilmicHejl          { "Filmic (Hejl)",          "filmic_hejl" };

    // TODO choose a single one after testing different presets
    // constexpr const ToneMapOperator FilmicUncharted     { "Filmic (Uncharted)",     "filmic_uncharted" };
    constexpr const ToneMapOperator FilmicUncharted1    { "Filmic (Uncharted) 1",   "filmic_uncharted1" };
    constexpr const ToneMapOperator FilmicUncharted2    { "Filmic (Uncharted) 2",   "filmic_uncharted2" };
    constexpr const ToneMapOperator FilmicUncharted3    { "Filmic (Uncharted) 3",   "filmic_uncharted3" };

    constexpr const ToneMapOperator Piecewise           { "Piecewise",              "piecewise" };
    constexpr const ToneMapOperator Reinhard            { "Reinhard",               "reinhard" };
    constexpr const ToneMapOperator ReinhardExtended    { "Reinhard (Extended)",    "reinhard_extended" };
    constexpr const ToneMapOperator DebugToneMap        { "Debug",                  "debug" };

    //@Todo add new TMOs
    #define TONE_MAP_OPERATOR_ARRAY {   \
        AcesNarkowicz.id,               \
        AcesUnreal.id,                  \
        FilmicHejl.id,                  \
        FilmicUncharted1.id,             \
        FilmicUncharted2.id,             \
        FilmicUncharted3.id,             \
        Piecewise.id,                   \
        Reinhard.id,                    \
        ReinhardExtended.id,            \
        DebugToneMap.id,                \
    }

    #define INSERT_TONE_MAP_OPERATOR(tmo) insert(tmo.label, tmo.id)

    //@Todo add new TMOs
    #define TONE_MAP_OPERATOR_DICTIONARY Dictionary()   \
        .INSERT_TONE_MAP_OPERATOR(AcesNarkowicz)        \
        .INSERT_TONE_MAP_OPERATOR(AcesUnreal)           \
        .INSERT_TONE_MAP_OPERATOR(FilmicHejl)           \
        .INSERT_TONE_MAP_OPERATOR(FilmicUncharted1)      \
        .INSERT_TONE_MAP_OPERATOR(FilmicUncharted2)      \
        .INSERT_TONE_MAP_OPERATOR(FilmicUncharted3)      \
        .INSERT_TONE_MAP_OPERATOR(Piecewise)            \
        .INSERT_TONE_MAP_OPERATOR(Reinhard)             \
        .INSERT_TONE_MAP_OPERATOR(ReinhardExtended)     \
        .INSERT_TONE_MAP_OPERATOR(DebugToneMap)

    //
    // Tone map post-processing stage.
    //

    const char* Model = "tone_map_post_processing_stage";

    constexpr const char* DeafaultToneMapOperatorId = AcesNarkowicz.id;

    class ToneMapPostProcessingStage
      : public PostProcessingStage
    {
      public:
        ToneMapPostProcessingStage(
            const char*             name,
            const ParamArray&       params)
          : PostProcessingStage(name, params)
          , m_tone_map(nullptr)
        {
        }

        void release() override
        {
            if (m_tone_map != nullptr)
                delete m_tone_map;

            delete this;
        }

        const char* get_model() const override
        {
            return Model;
        }

        bool on_frame_begin(
            const Project&          project,
            const BaseGroup*        parent,
            OnFrameBeginRecorder&   recorder,
            IAbortSwitch*           abort_switch) override
        {
            const OnFrameBeginMessageContext context("post-processing stage", this);

            const std::string tone_map_operator =
                m_params.get_optional<std::string>(
                    "tone_map_operator",
                    DeafaultToneMapOperatorId,
                    TONE_MAP_OPERATOR_ARRAY,
                    context);

            m_clip_values = m_params.get_optional("clip_values", true, context);

            //@Todo add new TMOs

            #define GET_OPT(name, default_value) m_params.get_optional(name, default_value, context)

            // Initialize the tone map applier.
            if (tone_map_operator == AcesNarkowicz.id)
            {
                m_tone_map = new AcesNarkowiczApplier();
            }
            else if (tone_map_operator == AcesUnreal.id)
            {
                m_tone_map = new AcesUnrealApplier();
            }
            else if (tone_map_operator == FilmicHejl.id)
            {
                m_tone_map = new FilmicHejlApplier();
            }
#if 0
            else if (tone_map_operator == FilmicUncharted.id)
            {
                const float A =
                    m_params.get_optional("filmic_uncharted_A", FilmicUnchartedApplier::DefaultA, context);
                const float B =
                    m_params.get_optional("filmic_uncharted_B", FilmicUnchartedApplier::DefaultB, context);
                const float C =
                    m_params.get_optional("filmic_uncharted_C", FilmicUnchartedApplier::DefaultC, context);
                const float D =
                    m_params.get_optional("filmic_uncharted_D", FilmicUnchartedApplier::DefaultD, context);
                const float E =
                    m_params.get_optional("filmic_uncharted_E", FilmicUnchartedApplier::DefaultE, context);
                const float F =
                    m_params.get_optional("filmic_uncharted_F", FilmicUnchartedApplier::DefaultF, context);
                const float W =
                    m_params.get_optional("filmic_uncharted_W", FilmicUnchartedApplier::DefaultW, context);
                const float exposure_bias =
                    m_params.get_optional("filmic_uncharted_exposure_bias", FilmicUnchartedApplier::DefaultExposureBias, context);

                m_tone_map = new FilmicUnchartedApplier(A, B, C, D, E, F, W, exposure_bias);
            }
#else
            else if (tone_map_operator == FilmicUncharted1.id)
            {
                m_tone_map =
                    new FilmicUnchartedApplier(
                        GET_OPT("filmic_uncharted1_A", FilmicUnchartedApplier::DefaultA),
                        GET_OPT("filmic_uncharted1_B", FilmicUnchartedApplier::DefaultB),
                        GET_OPT("filmic_uncharted1_C", FilmicUnchartedApplier::DefaultC),
                        GET_OPT("filmic_uncharted1_D", FilmicUnchartedApplier::DefaultD),
                        GET_OPT("filmic_uncharted1_E", FilmicUnchartedApplier::DefaultE),
                        GET_OPT("filmic_uncharted1_F", FilmicUnchartedApplier::DefaultF),
                        GET_OPT("filmic_uncharted1_W", FilmicUnchartedApplier::DefaultW),
                        GET_OPT("filmic_uncharted1_exposure_bias", FilmicUnchartedApplier::DefaultExposureBias));
            }
            else if (tone_map_operator == FilmicUncharted2.id)
            {
                m_tone_map =
                    new FilmicUnchartedApplier(
                        GET_OPT("filmic_uncharted2_A", FilmicUnchartedApplier::DefaultA),
                        GET_OPT("filmic_uncharted2_B", FilmicUnchartedApplier::DefaultB),
                        GET_OPT("filmic_uncharted2_C", FilmicUnchartedApplier::DefaultC),
                        GET_OPT("filmic_uncharted2_D", FilmicUnchartedApplier::DefaultD),
                        GET_OPT("filmic_uncharted2_E", FilmicUnchartedApplier::DefaultE),
                        GET_OPT("filmic_uncharted2_F", FilmicUnchartedApplier::DefaultF),
                        GET_OPT("filmic_uncharted2_W", FilmicUnchartedApplier::DefaultW),
                        GET_OPT("filmic_uncharted2_exposure_bias", FilmicUnchartedApplier::DefaultExposureBias));
            }
            else if (tone_map_operator == FilmicUncharted3.id)
            {
                m_tone_map =
                    new FilmicUnchartedApplier(
                        GET_OPT("filmic_uncharted3_A", FilmicUnchartedApplier::DefaultA),
                        GET_OPT("filmic_uncharted3_B", FilmicUnchartedApplier::DefaultB),
                        GET_OPT("filmic_uncharted3_C", FilmicUnchartedApplier::DefaultC),
                        GET_OPT("filmic_uncharted3_D", FilmicUnchartedApplier::DefaultD),
                        GET_OPT("filmic_uncharted3_E", FilmicUnchartedApplier::DefaultE),
                        GET_OPT("filmic_uncharted3_F", FilmicUnchartedApplier::DefaultF),
                        GET_OPT("filmic_uncharted3_W", FilmicUnchartedApplier::DefaultW),
                        GET_OPT("filmic_uncharted3_exposure_bias", FilmicUnchartedApplier::DefaultExposureBias));
            }
#endif
            else if (tone_map_operator == Piecewise.id)
            {
                const float toe_strength =
                    m_params.get_optional("piecewise_toe_strength", PiecewiseApplier::DefaultToeStrength, context);
                const float toe_length =
                    m_params.get_optional("piecewise_toe_length", PiecewiseApplier::DefaultToeLength, context);
                const float shoulder_strength =
                    m_params.get_optional("piecewise_shoulder_strength", PiecewiseApplier::DefaultShoulderStrength, context);
                const float shoulder_length =
                    m_params.get_optional("piecewise_shoulder_length", PiecewiseApplier::DefaultShoulderLength, context);
                const float shoulder_angle =
                    m_params.get_optional("piecewise_shoulder_angle", PiecewiseApplier::DefaultShoulderAngle, context);

                m_tone_map = new PiecewiseApplier(toe_strength, toe_length, shoulder_strength, shoulder_length, shoulder_angle);
            }
            else if (tone_map_operator == Reinhard.id)
            {
                const bool use_luminance =
                    m_params.get_optional("reinhard_use_luminance", ReinhardApplier::DefaultUseLuminance, context);

                m_tone_map = new ReinhardApplier(use_luminance);
            }
            else if (tone_map_operator == ReinhardExtended.id)
            {
                const float max_white =
                    m_params.get_optional("reinhard_extended_max_white", ReinhardExtendedApplier::DefaultMaxWhite, context);
                const bool use_luminance =
                    m_params.get_optional("reinhard_extended_use_luminance", ReinhardExtendedApplier::DefaultUseLuminance, context);

                m_tone_map = new ReinhardExtendedApplier(max_white, use_luminance);
            }
            else if (tone_map_operator == DebugToneMap.id)
            {
                m_tone_map = new DebugToneMapApplier();
            }
            else
            {
                // FIXME we shouldn't reach here.. but what if we do?
                // m_tone_map = nullptr;
                assert(false);
            }

            return true;
        }

        void execute(Frame& frame, const std::size_t thread_count) const override
        {
            const CanvasProperties& props = frame.image().properties();
            Image& image = frame.image();

            // Apply the selected tone mapping operator to each image tile, in parallel.
            m_tone_map->apply_on_tiles(image, thread_count);

            if (m_clip_values)
            {
                for (std::size_t y = 0; y < props.m_canvas_height; ++y)
                {
                    for (std::size_t x = 0; x < props.m_canvas_width; ++x)
                    {
                        Color3f color; // RGB
                        image.get_pixel(x, y, color);
                        image.set_pixel(x, y, saturate(color));
                    }
                }
            }
        }

      private:
        ToneMapApplier*     m_tone_map;
        bool                m_clip_values;
    };
}


//
// ToneMapPostProcessingStageFactory class implementation.
//

void ToneMapPostProcessingStageFactory::release()
{
    delete this;
}

const char* ToneMapPostProcessingStageFactory::get_model() const
{
    return Model;
}

Dictionary ToneMapPostProcessingStageFactory::get_model_metadata() const
{
    return
        Dictionary()
            .insert("name", Model)
            .insert("label", "Tone Map");
}

namespace
{
    inline void add_numeric_param_metadata(
        DictionaryArray&    metadata,
        const char*         name,
        const char*         label,
        const char*         min_value,
        const char*         min_type,
        const char*         max_value,
        const char*         max_type,
        const char*         default_value,
        const char*         tone_map_operator_id)
    {
        metadata.push_back(
            Dictionary()
                .insert("name", name)
                .insert("label", label)
                .insert("type", "numeric")
                .insert("min",
                        Dictionary()
                            .insert("value", min_value)
                            .insert("type", min_type))
                .insert("max",
                        Dictionary()
                            .insert("value", max_value)
                            .insert("type", max_type))
                .insert("use", "optional")
                .insert("default", default_value)
                .insert("visible_if",
                        Dictionary()
                            .insert("tone_map_operator", tone_map_operator_id)));
    }

    inline void add_boolean_param_metadata(
        DictionaryArray&    metadata,
        const char*         name,
        const char*         label,
        const char*         default_value,
        const char*         tone_map_operator_id)
    {
        metadata.push_back(
            Dictionary()
                .insert("name", name)
                .insert("label", label)
                .insert("type", "boolean")
                .insert("use", "optional")
                .insert("default", default_value)
                .insert("visible_if",
                        Dictionary()
                            .insert("tone_map_operator", tone_map_operator_id)));
    }
}

DictionaryArray ToneMapPostProcessingStageFactory::get_input_metadata() const
{
    DictionaryArray metadata;

    add_common_input_metadata(metadata);

    //@Todo remove (?)
    metadata.push_back(
        Dictionary()
            .insert("name", "clip_values")
            .insert("label", "Clip values")
            .insert("type", "boolean")
            .insert("use", "optional")
            .insert("default", "true"));

    metadata.push_back(
        Dictionary()
            .insert("name", "tone_map_operator")
            .insert("label", "Operator")
            .insert("type", "enumeration")
            .insert("items", TONE_MAP_OPERATOR_DICTIONARY)
            .insert("use", "required")
            .insert("default", DeafaultToneMapOperatorId)
            .insert("on_change", "rebuild_form"));

    //@Todo add new TMO params

    // ACES (Narkowicz)
    {
        // No parameters.
    }

    // ACES (Unreal)
    {
        // No parameters.
    }

    // Filmic (Hejl)
    {
        // No parameters.
    }

#if 0
    // Filmic (Uncharted)
    {
        add_numeric_param_metadata(
            metadata,
            "filmic_uncharted_A",
            "A (shoulder strength)",
            "0.0", "hard",              // min
            "1.0", "hard",              // max
            "0.22",                     // FilmicUnchartedApplier::DefaultA
            FilmicUncharted.id);

        add_numeric_param_metadata(
            metadata,
            "filmic_uncharted_B",
            "B (linear strength)",
            "0.0", "hard",              // min
            "1.0", "hard",              // max
            "0.30",                     // FilmicUnchartedApplier::DefaultB
            FilmicUncharted.id);

        add_numeric_param_metadata(
            metadata,
            "filmic_uncharted_C",
            "C (linear angle)",
            "0.0", "hard",              // min
            "1.0", "hard",              // max
            "0.10",                     // FilmicUnchartedApplier::DefaultC
            FilmicUncharted.id);

        add_numeric_param_metadata(
            metadata,
            "filmic_uncharted_D",
            "D (toe strength)",
            "0.0", "hard",              // min
            "1.0", "hard",              // max
            "0.20",                     // FilmicUnchartedApplier::DefaultD
            FilmicUncharted.id);

        add_numeric_param_metadata(
            metadata,
            "filmic_uncharted_E",
            "E (toe numerator)",
            "0.0", "hard",              // min
            "1.0", "hard",              // max
            "0.01",                     // FilmicUnchartedApplier::DefaultE
            FilmicUncharted.id);

        add_numeric_param_metadata(
            metadata,
            "filmic_uncharted_F",
            "F (toe denominator)",
            "0.0", "hard",              // min
            "1.0", "hard",              // max
            "0.30",                     // FilmicUnchartedApplier::DefaultF
            FilmicUncharted.id);

        add_numeric_param_metadata(
            metadata,
            "filmic_uncharted_W",
            "Linear white point",
            "0.0", "hard",              // min
            "1.0", "hard",              // max
            "11.2",                     // FilmicUnchartedApplier::DefaultW
            FilmicUncharted.id);

        add_numeric_param_metadata(
            metadata,
            "filmic_uncharted_exposure_bias",
            "Exposure bias",
            "0.0", "hard",              // min
            "10.0", "hard",             // max
            "2.0",                      // FilmicUnchartedApplier::DefaultExposureBias
            FilmicUncharted.id);
    }
#else
    // Filmic (Uncharted) 1
    {
        #define ADD_PARAM_1(name, label, max_value, default_value) \
            add_numeric_param_metadata(metadata, name, label, "0.0", "hard", max_value, "hard", default_value, FilmicUncharted1.id);

        // Original presentation values:
        ADD_PARAM_1("filmic_uncharted1_A", "A (shoulder strength)", "1.0", "0.22");
        ADD_PARAM_1("filmic_uncharted1_B", "B (linear strength)", "1.0", "0.30");
        ADD_PARAM_1("filmic_uncharted1_C", "C (linear angle)", "1.0", "0.10");
        ADD_PARAM_1("filmic_uncharted1_D", "D (toe strength)", "1.0", "0.20");
        ADD_PARAM_1("filmic_uncharted1_E", "E (toe numerator)", "1.0", "0.01");
        ADD_PARAM_1("filmic_uncharted1_F", "F (toe denominator)", "1.0", "0.30");
        ADD_PARAM_1("filmic_uncharted1_W", "Linear white point", "1.0", "11.2");
        ADD_PARAM_1("filmic_uncharted1_exposure_bias", "Exposure bias", "10.0", "2.0");
    }

    // Filmic (Uncharted) 2
    {
        #define ADD_PARAM_2(name, label, max_value, default_value) \
            add_numeric_param_metadata(metadata, name, label, "0.0", "hard", max_value, "hard", default_value, FilmicUncharted2.id);

        // Updated values (on follow up blog):
        ADD_PARAM_2("filmic_uncharted2_A", "A (shoulder strength)", "1.0", "0.15");
        ADD_PARAM_2("filmic_uncharted2_B", "B (linear strength)", "1.0", "0.50");
        ADD_PARAM_2("filmic_uncharted2_C", "C (linear angle)", "1.0", "0.10");
        ADD_PARAM_2("filmic_uncharted2_D", "D (toe strength)", "1.0", "0.20");
        ADD_PARAM_2("filmic_uncharted2_E", "E (toe numerator)", "1.0", "0.02");
        ADD_PARAM_2("filmic_uncharted2_F", "F (toe denominator)", "1.0", "0.30");
        ADD_PARAM_2("filmic_uncharted2_W", "Linear white point", "1.0", "11.2");
        ADD_PARAM_2("filmic_uncharted2_exposure_bias", "Exposure bias", "10.0", "2.0");
    }

    // Filmic (Uncharted) 3
    {
        #define ADD_PARAM_3(name, label, max_value, default_value) \
            add_numeric_param_metadata(metadata, name, label, "0.0", "hard", max_value, "hard", default_value, FilmicUncharted3.id);

        // Values used in MJP's BakingLab
        ADD_PARAM_3("filmic_uncharted3_A", "A (shoulder strength)", "10.0", "4.0");
        ADD_PARAM_3("filmic_uncharted3_B", "B (linear strength)", "10.0", "5.0");
        ADD_PARAM_3("filmic_uncharted3_C", "C (linear angle)", "1.0", "0.12");
        ADD_PARAM_3("filmic_uncharted3_D", "D (toe strength)", "20.0", "13.0");
        ADD_PARAM_3("filmic_uncharted3_E", "E (toe numerator)", "1.0", "0.01");
        ADD_PARAM_3("filmic_uncharted3_F", "F (toe denominator)", "1.0", "0.30");
        ADD_PARAM_3("filmic_uncharted3_W", "Linear white point", "1.0", "6.0");
        ADD_PARAM_3("filmic_uncharted3_exposure_bias", "Exposure bias", "10.0", "1.0");
    }
#endif

    // Piecewise
    {
        add_numeric_param_metadata(
            metadata,
            "piecewise_toe_strength",
            "Toe Strength",
            "0.0", "soft",              // min
            "1.0", "soft",              // max
            "0.0",                      // PiecewiseApplier::DefaultToeStrength
            Piecewise.id);

        add_numeric_param_metadata(
            metadata,
            "piecewise_toe_length",
            "Toe Length",
            "0.0", "soft",              // min
            "1.0", "soft",              // max
            "0.5",                      // PiecewiseApplier::DefaultToeLength
            Piecewise.id);

        add_numeric_param_metadata(
            metadata,
            "piecewise_shoulder_strength",
            "Shoulder Strength",
            "0.0", "soft",              // min
            "1.0", "soft",              // max
            "0.0",                      // PiecewiseApplier::DefaultShoulderStrength
            Piecewise.id);

        add_numeric_param_metadata(
            metadata,
            "piecewise_shoulder_length",
            "Shoulder Length (F-stops)",

            // FIXME this is expressed in F-stops instead of
            // as a ratio, thus, what should be its min/max?
            "0.00001", "hard",          // min
            "10.0", "soft",             // max

            "0.5",                      // PiecewiseApplier::DefaultShoulderLength
            Piecewise.id);

        add_numeric_param_metadata(
            metadata,
            "piecewise_shoulder_angle",
            "Shoulder Angle",
            "0.0", "soft",              // min
            "1.0", "soft",              // max
            "0.0",                      // PiecewiseApplier::DefaultShoulderAngle
            Piecewise.id);
    }

    // Reinhard
    {
        add_boolean_param_metadata(
            metadata,
            "reinhard_use_luminance",
            "Use Luminance",
            "true",                     // ReinhardApplier::DefaultUseLuminance
            Reinhard.id);
    }

    // Reinhard (Extended)
    {
        add_numeric_param_metadata(
            metadata,
            "reinhard_extended_max_white",
            "Lmax",

            // FIXME min/max luminance values can only
            // be accurately computed at run-time.. :(
            "0.0", "hard",              // min
            "10000.0", "soft",          // max

            "1.0",                      // ReinhardExtendedApplier::DefaultMaxWhite
            ReinhardExtended.id);

        add_boolean_param_metadata(
            metadata,
            "reinhard_extended_use_luminance",
            "Use Luminance",
            "true",                     // ReinhardExtendedApplier::DefaultUseLuminance
            ReinhardExtended.id);
    }

    return metadata;
}

auto_release_ptr<PostProcessingStage> ToneMapPostProcessingStageFactory::create(
    const char*         name,
    const ParamArray&   params) const
{
    return auto_release_ptr<PostProcessingStage>(
        new ToneMapPostProcessingStage(name, params));
}

}   // namespace renderer
