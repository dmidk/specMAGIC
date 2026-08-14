#pragma once

#include <algorithm>   // std::copy, std::fill
#include <stdexcept>   // std::runtime_error, std::out_of_range
#include <utility>     // std::swap
#include <cstdint>     // std::int64_t
#include "types.hpp"   

struct Matrix {
    
    public:

        // Construct a rows x cols matrix, indexed [0..rows-1][0..cols-1]
        Matrix(int rows, int cols)
            : rows_(rows), cols_(cols), data_(nullptr)
        {
            if (rows_ < 0 || cols_ < 0) {
                throw std::runtime_error("Matrix: negative dimension");
            }

            // Avoid overflow in rows*cols 
            const std::int64_t count64 = static_cast<std::int64_t>(rows_) * static_cast<std::int64_t>(cols_);
            if (count64 < 0) {
                throw std::runtime_error("Matrix: size overflow");
            }

            const int count = static_cast<int>(count64);
            data_ = (count > 0) ? new MAGIC_INT[count]{} : nullptr; // zero-initialise
        }

        // Default ctor
        // Probably don't need this, just for consistency
        Matrix() noexcept = default;

        // Dtor
        // Since we use dynamic allocs
        ~Matrix() {
            delete[] data_;
        }

        // --------------------------------------------
        // All of the operators in this block are probably totally overkill to define
        // But doing so avoids weird errors later on

        Matrix(const Matrix&) = delete;
        Matrix& operator=(const Matrix&) = delete;  

        // Move constructor
        Matrix(Matrix&& other) noexcept
            : rows_(other.rows_), cols_(other.cols_), data_(other.data_)
        {
            other.rows_ = other.cols_ = 0;
            other.data_ = nullptr;
        }

        // Move assignment
        Matrix& operator=(Matrix&& other) noexcept {
            if (this != &other) {
                delete[] data_;
                rows_ = other.rows_;
                cols_ = other.cols_;
                data_ = other.data_;

                other.rows_ = other.cols_ = 0;
                other.data_ = nullptr;
            }
            return *this;
        }

        void swap(Matrix& other) noexcept {
            std::swap(rows_, other.rows_);
            std::swap(cols_, other.cols_);
            std::swap(data_, other.data_);
        }
        // ---------------------------------------------------------

        // Below we define both const and non const overloads
        // Probably totally overkill
        // But avoids compiler nonsense
        
        // Unchecked access. Assumes 0 <= i < rows and 0 <= j < cols.
        // Prefer not to use these.
        MAGIC_INT& operator()(int i, int j) noexcept {
            return data_[i * cols_ + j];
        }
        const MAGIC_INT& operator()(int i, int j) const noexcept {
            return data_[i * cols_ + j];
        }

        // Checked access (throws std::out_of_range)
        MAGIC_INT& at(int i, int j) {
            if (i < 0 || i >= rows_ || j < 0 || j >= cols_) {
                throw std::out_of_range("Matrix::at: index out of range");
            }
            return data_[i * cols_ + j];
        }
        const MAGIC_INT& at(int i, int j) const {
            if (i < 0 || i >= rows_ || j < 0 || j >= cols_) {
                throw std::out_of_range("Matrix::at: index out of range");
            }
            return data_[i * cols_ + j];
        }

        // Getters
        MAGIC_INT* data() noexcept { return data_; }
        const MAGIC_INT* data() const noexcept { return data_; }

        int rows() const noexcept { return rows_; }
        int cols() const noexcept { return cols_; }


    private:
        int rows_{0}, cols_{0};
        MAGIC_INT* data_{nullptr};
};