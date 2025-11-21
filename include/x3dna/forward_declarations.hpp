/**
 * @file forward_declarations.hpp
 * @brief Forward declarations for X3DNA classes
 */

#pragma once

namespace x3dna {

// Core classes
namespace core {
    class Atom;
    class Residue;
    class Chain;
    class Structure;
    class BasePair;
    class ReferenceFrame;
    struct BasePairStepParameters;
    struct HelicalParameters;
}

// Geometry classes
namespace geometry {
    class Vector3D;
    class Matrix3D;
    class Quaternion;
    class LeastSquaresFitter;
}

// I/O classes
namespace io {
    class PdbParser;
    class PdbWriter;
    class JsonReader;
    class JsonWriter;
    class InputFileParser;
}

// Algorithm classes
namespace algorithms {
    class BaseFrameCalculator;
    class BasePairFinder;
    class ParameterCalculator;
    class HelixDetector;
    class HydrogenBondValidator;
    class StandardBaseTemplates;
    class RingAtomMatcher;
}

// Protocol classes
namespace protocols {
    class ProtocolBase;
    class FindPairProtocol;
    class AnalyzeProtocol;
}

// Configuration
namespace config {
    class ConfigManager;
    struct ParameterThresholds;
}

} // namespace x3dna

